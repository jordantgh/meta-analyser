from sqlalchemy import create_engine, Column, String, text
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
from uuid import uuid4
import pandas as pd

Base = declarative_base()


class ProcessedTableDBEntry(Base):
    __tablename__ = 'processed_tables'
    table_id = Column(String, primary_key=True)
    original_file_id = Column(String)
    sql_table_name = Column(String)


class PostPruningTableDBEntry(Base):
    __tablename__ = 'post_pruning_tables'
    table_id = Column(String, primary_key=True)
    original_file_id = Column(String)
    sql_table_name = Column(String)


class TableDBManager:
    def __init__(
        self,
        processed_db_url=f"sqlite:///processed_tables-{str(uuid4())}.db",
        post_pruning_db_url=f"sqlite:///post_pruning_tables-{str(uuid4())}.db"
    ):
        self.processed_db_url = processed_db_url
        self.post_pruning_db_url = post_pruning_db_url
        
        self.processed_engine = create_engine(processed_db_url)
        self.post_pruning_engine = create_engine(post_pruning_db_url)

        Base.metadata.create_all(self.processed_engine)
        Base.metadata.create_all(self.post_pruning_engine)

        self.ProcessedSession = sessionmaker(bind=self.processed_engine)
        self.PostPruningSession = sessionmaker(bind=self.post_pruning_engine)

    def _get_engine_and_session(self, table_class):
        if table_class == ProcessedTableDBEntry:
            return self.processed_engine, self.ProcessedSession
        elif table_class == PostPruningTableDBEntry:
            return self.post_pruning_engine, self.PostPruningSession
        else:
            raise ValueError("Invalid table_class")

    def save_table(self, table_class, table_id, original_file_id, df):
        engine, Session = self._get_engine_and_session(table_class)
        with Session() as session:
            sql_table_name = f"{table_class.__tablename__}_{table_id}"
            df.to_sql(sql_table_name, engine, index=False)
            new_table = table_class(
                table_id=table_id,
                original_file_id=str(original_file_id),
                sql_table_name=sql_table_name
            )
            session.add(new_table)
            session.commit()

    def update_table(self, table_class, table_id, df):
        engine, Session = self._get_engine_and_session(table_class)
        with Session() as session:
            existing_table = session.query(table_class).filter_by(table_id=table_id).first()
            if existing_table:
                df.to_sql(existing_table.sql_table_name, engine, if_exists='replace', index=False)
                session.commit()

    def get_processed_table_data(self, table_id):
        return self.get_table_data(ProcessedTableDBEntry, table_id)

    def get_post_pruning_table_data(self, table_id):
        return self.get_table_data(PostPruningTableDBEntry, table_id)

    def get_table_data(self, table_class, table_id):
        engine, Session = self._get_engine_and_session(table_class)
        with Session() as session:
            table_entry = session.query(table_class).filter_by(table_id=table_id).first()
            if table_entry:
                df = pd.read_sql_table(table_entry.sql_table_name, engine)
                return df.reset_index(drop=True)
            return None

    def get_table_object(self, table_class, table_id):
        _, Session = self._get_engine_and_session(table_class)
        with Session() as session:
            table = session.query(table_class).filter_by(table_id=table_id).first()
            return table

    def delete_table(self, table_class, table_id):
        engine, Session = self._get_engine_and_session(table_class)
        with Session() as session:
            table_entry = session.query(table_class).filter_by(table_id=table_id).first()
            if table_entry:
                with engine.connect() as conn:
                    conn.execute(text(f"DROP TABLE IF EXISTS \"{table_entry.sql_table_name}\""))
                session.delete(table_entry)
                session.commit()

    def reset(self):
        for engine, Session in [
            (self.processed_engine, self.ProcessedSession),
            (self.post_pruning_engine, self.PostPruningSession)]:
            with Session() as session:
                for table_class in [ProcessedTableDBEntry, PostPruningTableDBEntry]:
                    table_entries = session.query(table_class).all()
                    for entry in table_entries:
                        with engine.connect() as conn:
                            conn.execute(text(f"DROP TABLE IF EXISTS \"{entry.sql_table_name}\""))
                    session.query(table_class).delete()
                session.commit()

def processed_df_to_db(db_manager, table_id, original_file_id, df):
    db_manager.save_table(ProcessedTableDBEntry, table_id, original_file_id, df)
