import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from API import *
from project import *
import time

base_dir = os.environ.get("BASE_DIR")
NPB_dir = f"{base_dir}/NPB3.0-omp-C"
file_dir = f'{base_dir}/NPB3.0-omp-C/MG/mg_#_omp.c'


questions = ["你好", "今天天气怎么样", "Python有什么优点","1+1=?","写一个冒泡排序"]

start = time.time()
for q in questions:
    print(ask_gpt_question(q))  # 假设每次1秒，10个问题需要10秒
print(f"同步耗时: {time.time() - start}")

# 异步：耗时约等于最慢的那个请求
async def test_async_batch():
    start = time.time()
    answers = await batch_ask_questions(questions)
    print(f"异步耗时: {time.time() - start}")
    print(answers)
    return answers

answers = asyncio.run(test_async_batch())

