from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker



mysql_psw = os.environ["MYSQL_PWD"]
sql_alchemy_database_url = "sqlite:///./data.db"
engine  = create_engine(sql_alchemy_database_url, connect_args={"check_same_thread": False}, echo=True)
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)

Base = declarative_base()

def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()