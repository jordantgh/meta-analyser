from typing import TYPE_CHECKING, Optional
if TYPE_CHECKING:
    from uuid import UUID
    from pandas import DataFrame
    from utils.constants import PageIdentity

from sqlalchemy import create_engine, Column, String, text
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.dialects.sqlite import JSON
from sqlalchemy.orm import sessionmaker
from uuid import uuid4
import pandas as pd
import shutil
import os

from utils.constants import PageIdentity

Base = declarative_base()


class TableDBEntry(Base):
    __abstract__ = True
    pmc_id = Column(String)
    paper_title = Column(String)
    paper_url = Column(String)
    supp_url = Column(String)
    table_id = Column(String, primary_key=True)
    tags = Column(JSON)


class ProcessedTableDBEntry(TableDBEntry):
    __tablename__ = 'processed_tables'


class PostPruningTableDBEntry(TableDBEntry):
    __tablename__ = 'post_pruning_tables'


class SortedListTableDBEntry(TableDBEntry):
    __tablename__ = 'sorted_list_tables'


class TableDBManager:
    def __init__(
        self,
        db_temp_path_root: 'str' = ".",
        db_perm_path_root: 'str' = ".",
        db_session_files: 'Optional[list[str]]' = None
    ):

        self.db_temp_path_root = db_temp_path_root
        self.db_perm_path_root = db_perm_path_root

        if db_session_files:
            self._create_temp_dbs(db_temp_path_root, db_session_files)
        else:
            self._create_temp_dbs(db_temp_path_root)

        self._establish_connections()

        Base.metadata.create_all(self.processed_engine)
        Base.metadata.create_all(self.post_pruning_engine)
        Base.metadata.create_all(self.sorted_list_engine)

    def _create_temp_dbs(
        self, db_path: 'str',
        db_session_files: 'Optional[list[str]]' = None
    ):
        id = str(uuid4())
        self.processed_db_temp_path = os.path.join(
            db_path, f"processed_db-{id}.db"
        )
        self.pruned_db_temp_path = os.path.join(db_path, f"pruned_db-{id}.db")
        self.sorted_list_db_temp_path = os.path.join(
            db_path, f"sorted_list_db-{id}.db")
        if db_session_files:
            shutil.copyfile(db_session_files[0], self.processed_db_temp_path)
            shutil.copyfile(db_session_files[1], self.pruned_db_temp_path)
            shutil.copyfile(db_session_files[2], self.sorted_list_db_temp_path)

    def _get_engine_and_session(self, table_class: 'TableDBEntry') -> 'tuple':
        if table_class == ProcessedTableDBEntry:
            return self.processed_engine, self.ProcessedSession
        elif table_class == PostPruningTableDBEntry:
            return self.post_pruning_engine, self.PostPruningSession
        elif table_class == SortedListTableDBEntry:
            return self.sorted_list_engine, self.SortedListSession
        else:
            raise ValueError("Invalid table_class")

    def _establish_connections(self):
        self.processed_engine = create_engine(
            f"sqlite:///{self.processed_db_temp_path}")
        self.post_pruning_engine = create_engine(
            f"sqlite:///{self.pruned_db_temp_path}")
        self.sorted_list_engine = create_engine(
            f"sqlite:///{self.sorted_list_db_temp_path}")
        self.ProcessedSession = sessionmaker(bind=self.processed_engine)
        self.PostPruningSession = sessionmaker(bind=self.post_pruning_engine)
        self.SortedListSession = sessionmaker(bind=self.sorted_list_engine)

    def _dispose_connections(self):
        self.processed_engine.dispose()
        self.post_pruning_engine.dispose()
        self.sorted_list_engine.dispose()

    def save_table(
        self,
        table_class: 'TableDBEntry',
        pmc_id: 'str',
        paper_title: 'str',
        paper_url: 'str',
        supp_url: 'str',
        df: 'DataFrame',
        table_id: 'str',
        tags: 'list[str]' = []
    ):
        engine, Session = self._get_engine_and_session(table_class)
        with Session() as session:
            df.to_sql(table_id, engine, index=False)
            new_table = table_class(
                pmc_id=pmc_id,
                paper_title=paper_title,
                paper_url=paper_url,
                supp_url=supp_url,
                table_id=table_id,
                tags=tags
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
                    existing_table.table_id,
                    engine,
                    if_exists='replace',
                    index=False
                )
                session.commit()

    def update_table_tags(self, table_id: 'str', tags: 'list[str]'):
        self._update_table_tags(ProcessedTableDBEntry, table_id, tags)
        self._update_table_tags(PostPruningTableDBEntry, table_id, tags)
        self._update_table_tags(SortedListTableDBEntry, table_id, tags)

    def _update_table_tags(
        self,
        table_class: 'TableDBEntry',
        table_id: 'str',
        tags: 'list[str]'
    ):
        _, Session = self._get_engine_and_session(table_class)
        with Session() as session:
            table_entry: 'TableDBEntry' = session.query(
                table_class).filter_by(table_id=table_id).first()
            if table_entry:
                table_entry.tags = tags
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
                df = pd.read_sql_table(table_entry.table_id, engine)
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
                        text(f"DROP TABLE IF EXISTS \"{table_entry.table_id}\""))
                session.delete(table_entry)
                session.commit()

    def save_dbs(
        self,
        filename: 'str',
        db_session_files: 'Optional[list[str]]' = None
    ) -> 'list[str]':

        if db_session_files:
            os.remove(db_session_files[0])
            os.remove(db_session_files[1])
            os.remove(db_session_files[2])

        basename = os.path.basename(filename)
        basename_no_ext = os.path.splitext(basename)[0]

        processed_target_path = os.path.join(
            self.db_perm_path_root,
            f"{basename_no_ext}-processed.db"
        )

        pruned_target_path = os.path.join(
            self.db_perm_path_root,
            f"{basename_no_ext}-pruned.db"
        )

        sorted_list_target_path = os.path.join(
            self.db_perm_path_root,
            f"{basename_no_ext}-sorted-list.db"
        )

        self._dispose_connections()

        shutil.copyfile(self.processed_db_temp_path, processed_target_path)
        shutil.copyfile(self.pruned_db_temp_path, pruned_target_path)
        shutil.copyfile(self.sorted_list_db_temp_path, sorted_list_target_path)

        self._establish_connections()

        return [processed_target_path, pruned_target_path, sorted_list_target_path]

    def reset(self):
        for engine, Session in [
            (self.processed_engine, self.ProcessedSession),
            (self.post_pruning_engine, self.PostPruningSession),
            (self.sorted_list_engine, self.SortedListSession)
        ]:
            with Session() as session:
                for table_class in [
                    ProcessedTableDBEntry,
                    PostPruningTableDBEntry,
                    SortedListTableDBEntry
                ]:
                    table_entries: 'list[TableDBEntry]' = session.query(
                        table_class
                    ).all()
                    for entry in table_entries:
                        with engine.connect() as conn:
                            conn.execute(
                                text(f"DROP TABLE IF EXISTS \"{entry.table_id}\""))
                    session.query(table_class).delete()
                session.commit()

    def delete_dbs(self):
        # Close all connections
        self.processed_engine.dispose()
        self.post_pruning_engine.dispose()
        self.sorted_list_engine.dispose()

        # Delete databases
        if os.path.exists(self.processed_db_temp_path):
            os.remove(self.processed_db_temp_path)

        if os.path.exists(self.pruned_db_temp_path):
            os.remove(self.pruned_db_temp_path)

        if os.path.exists(self.sorted_list_db_temp_path):
            os.remove(self.sorted_list_db_temp_path)


def processed_df_to_db(
    db_manager: 'TableDBManager',
    pmc_id: 'str',
    paper_title: 'str',
    paper_url: 'str',
    supp_url: 'str',
    df: 'DataFrame',
    table_id: 'str',
    tags: 'list[str]'
):
    db_manager.save_table(
        ProcessedTableDBEntry,
        pmc_id,
        paper_title,
        paper_url,
        supp_url,
        df,
        table_id,
        tags
    )
