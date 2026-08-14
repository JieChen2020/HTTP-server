from fastapi import FastAPI
from pydantic import BaseModel
from typing import List, Any, Dict
from bayesian_optimization import get_next_exps
from exp import run_exp
import uvicorn
import traceback

app = FastAPI()

class ExpData(BaseModel):
    condition: List
    outcomes: List
    q: int
    goal: List[Dict]
    num_of_init: int
    design_space: List[Dict]
    reaction: str

@app.post("/get-next-exps")
async def next_exps(data: ExpData):
    try:
        return get_next_exps(data.model_dump())
    except Exception as e:
        return {"error": str(e), "traceback": traceback.format_exc()}

@app.post("/run-exp")
async def exp(data: ExpData):
    try:
        return run_exp(data.model_dump())
    except Exception as e:
        return {"error": str(e), "traceback": traceback.format_exc()}

if __name__ == '__main__':
    uvicorn.run(app="http_api:app", host="10.97.40.181", port=123)
