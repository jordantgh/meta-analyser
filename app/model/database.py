from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from uuid import UUID
    from pandas import DataFrame
    from utils.constants import PageIdentity

from sqlalchemy import create_engine, Column, String, text
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
from uuid import uuid4
import pandas as pd
import shutil

from utils.constants import PageIdentity

Base = declarative_base()


class TableDBEntry(Base):
    __abstract__ = True
    table_id = Column(String, primary_key=True)
    original_file_id = Column(String)
    sql_table_name = Column(String)


class ProcessedTableDBEntry(TableDBEntry):
    __tablename__ = 'processed_tables'


class PostPruningTableDBEntry(TableDBEntry):
    __tablename__ = 'post_pruning_tables'


class TableDBManager:
    def __init__(
        self,
        db_path: 'str' = ".",
        processed_db_url: 'str' = None,
        post_pruning_db_url: 'str' = None
    ):
        self.db_path = db_path

        if processed_db_url and post_pruning_db_url:
            self.processed_db_url = processed_db_url
            self.post_pruning_db_url = post_pruning_db_url
        else:
            self._create_db_urls(db_path)

        self.processed_engine = create_engine(self.processed_db_url)
        self.post_pruning_engine = create_engine(self.post_pruning_db_url)

        Base.metadata.create_all(self.processed_engine)
        Base.metadata.create_all(self.post_pruning_engine)

        self.ProcessedSession = sessionmaker(bind=self.processed_engine)
        self.PostPruningSession = sessionmaker(bind=self.post_pruning_engine)

    def _create_db_urls(self, db_path: 'str'):
        id = str(uuid4())
        processed_db_url = (
            f"sqlite:///{db_path}/processed_tables-{id}.db"
        )

        post_pruning_db_url = (
            f"sqlite:///{db_path}/post_pruning_tables-{id}.db"
        )

        self.processed_db_url = processed_db_url
        self.post_pruning_db_url = post_pruning_db_url
    
    def _get_engine_and_session(self, table_class: 'TableDBEntry') -> 'tuple':
        if table_class == ProcessedTableDBEntry:
            return self.processed_engine, self.ProcessedSession
        elif table_class == PostPruningTableDBEntry:
            return self.post_pruning_engine, self.PostPruningSession
        else:
            raise ValueError("Invalid table_class")

    def save_table(
        self,
        table_class: 'TableDBEntry',
        table_id: 'str',
        original_file_id: 'UUID',
        df: 'DataFrame'
    ):
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

    def update_table(
        self,
        table_class: 'TableDBEntry',
        table_id: 'str',
        df: 'DataFrame'
    ):
        engine, Session = self._get_engine_and_session(table_class)
        with Session() as session:
            existing_table: 'TableDBEntry' = session.query(
                table_class
            ).filter_by(table_id=table_id).first()

            if existing_table:
                df.to_sql(
                    existing_table.sql_table_name,
                    engine,
                    if_exists='replace',
                    index=False
                )
                session.commit()

    def get_processed_table_data(
        self,
        table_id: 'str',
        context: 'PageIdentity'
    ) -> 'DataFrame':
        if context == PageIdentity.PARSED:
            return self.get_table_data(ProcessedTableDBEntry, table_id)
        else:
            return self.get_table_data(PostPruningTableDBEntry, table_id)

    def get_table_data(
        self,
        table_class: 'TableDBEntry',
        table_id: 'str'
    ) -> 'DataFrame':
        engine, Session = self._get_engine_and_session(table_class)
        with Session() as session:
            table_entry: 'TableDBEntry' = session.query(table_class).filter_by(
                table_id=table_id
            ).first()
            if table_entry:
                df = pd.read_sql_table(table_entry.sql_table_name, engine)
                return df.reset_index(drop=True)
            return None

    def get_table_object(
        self,
        table_class: 'TableDBEntry',
        table_id: 'str'
    ) -> 'TableDBEntry':
        _, Session = self._get_engine_and_session(table_class)
        with Session() as session:
            table = session.query(table_class).filter_by(
                table_id=table_id).first()
            return table

    def delete_table(
        self,
        table_class: 'TableDBEntry',
        table_id: 'str'
    ):
        engine, Session = self._get_engine_and_session(table_class)
        with Session() as session:
            table_entry: 'TableDBEntry' = session.query(table_class).filter_by(
                table_id=table_id).first()
            if table_entry:
                with engine.connect() as conn:
                    conn.execute(
                        text(f"DROP TABLE IF EXISTS \"{table_entry.sql_table_name}\""))
                session.delete(table_entry)
                session.commit()

    def copy_dbs(self) -> 'list[str]':
        new_id = str(uuid4())

        processed_copy = f"{self.db_path}/processed_tables-{new_id}.db"
        post_pruning_copy = f"{self.db_path}/post_pruning_tables-{new_id}.db"

        processed_src = self.processed_db_url.replace("sqlite:///", "")
        post_pruning_src = self.post_pruning_db_url.replace("sqlite:///", "")
        
        shutil.copyfile(processed_src, processed_copy)
        shutil.copyfile(post_pruning_src, post_pruning_copy)

        return [f"sqlite:///{p}" for p in [processed_copy, post_pruning_copy]]

    def reset(self):
        for engine, Session in [
            (self.processed_engine, self.ProcessedSession),
                (self.post_pruning_engine, self.PostPruningSession)]:
            with Session() as session:
                for table_class in [
                    ProcessedTableDBEntry, PostPruningTableDBEntry
                ]:
                    table_entries: 'list[TableDBEntry]' = session.query(
                        table_class
                    ).all()
                    for entry in table_entries:
                        with engine.connect() as conn:
                            conn.execute(
                                text(f"DROP TABLE IF EXISTS \"{entry.sql_table_name}\""))
                    session.query(table_class).delete()
                session.commit()


def processed_df_to_db(
    db_manager: 'TableDBManager',
    table_id: 'str',
    original_file_id: 'UUID',
    df: 'DataFrame'
):
    db_manager.save_table(ProcessedTableDBEntry, table_id, original_file_id, df)
