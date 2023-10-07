from sqlalchemy import create_engine, Column, String, BLOB
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
from uuid import uuid4
import pickle

Base = declarative_base()


class ProcessedTableDBEntry(Base):
    __tablename__ = 'processed_tables'

    table_id = Column(String, primary_key=True)
    original_file_id = Column(String)
    table_data = Column(BLOB)


class PostPruningTableDBEntry(Base):
    __tablename__ = 'post_pruning_tables'

    table_id = Column(String, primary_key=True)
    original_file_id = Column(String)
    table_data = Column(BLOB)


class TableDBManager:
    def __init__(self, db_url=f"sqlite:///tables-{str(uuid4())}.db"):
        self.engine = create_engine(db_url)
        Base.metadata.create_all(self.engine)
        self.Session = sessionmaker(bind=self.engine)

    def save_table(self, table_class, table_id, original_file_id, table_data):
        with self.Session() as session:
            new_table = table_class(
                table_id=table_id,
                original_file_id=str(original_file_id),
                table_data=table_data
            )
            session.add(new_table)
            session.commit()

    def update_table(self, table_class, table_id, table_data):
        with self.Session() as session:
            existing_table = session.query(
                table_class).filter_by(table_id=table_id).first()
            if existing_table:
                existing_table.table_data = table_data
                session.commit()

    def get_processed_table_data(self, table_id):
        return self.get_table_data(ProcessedTableDBEntry, table_id)

    def get_post_pruning_table_data(self, table_id):
        return self.get_table_data(PostPruningTableDBEntry, table_id)

    def get_table_data(self, table_class, table_id):
        with self.Session() as session:
            table = session.query(table_class).filter_by(
                table_id=table_id).first()
            if table:
                table_data = pickle.loads(table.table_data)
                table_data.reset_index(drop=True, inplace=True)
                return table_data
            return None

    def get_table_object(self, table_class, table_id):
        with self.Session() as session:
            table = session.query(table_class).filter_by(
                table_id=table_id).first()
            return table

    def delete_table(self, table_class, table_id):
        with self.Session() as session:
            table = session.query(table_class).filter_by(
                table_id=table_id).first()
            if table:
                session.delete(table)
                session.commit()

    def reset(self):
        with self.Session() as session:
            session.query(ProcessedTableDBEntry).delete()
            session.query(PostPruningTableDBEntry).delete()
            session.commit()


def processed_df_to_db(db_manager, table_id, original_file_id, df):
    serialized_df = pickle.dumps(df)
    db_manager.save_table(ProcessedTableDBEntry, table_id,
                          original_file_id, serialized_df)
