from fastapi import FastAPI, Depends, HTTMException
from src.db import engine, get_db
from sqlalchemy.orm import Session
from src.create_database import *
from typing import List
from fastapi import Path


@app.get("/genes/{chromosome}:{start}-{end}/", response_model=List[Gencode_genes])