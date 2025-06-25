import openai
import time
import re
import os
import subprocess
import asyncio
from openai import AsyncOpenAI
base_dir = os.environ.get("BASE_DIR")
base_url = 'https://open.xiaojingai.com/v1/'
api_key = ''
# openai.base_url = 'https://api.siliconflow.cn/v1/'
# openai.api_key = ''

openai.base_url = base_url
openai.api_key = api_key
client = AsyncOpenAI(
    base_url = base_url,
    api_key = api_key
)


def ask_gpt_question(prompt, max_retries=5, retry_delay=5):
    """
    最简单的API调用, print(ask_gpt_question("hello")
    """
    retries = 0
    while retries < max_retries:
        try:
            response = openai.chat.completions.create(
                model='o4-mini-high-all',  
                messages=[
                    {"role": "user", "content": prompt}
                ],
            )
            return response.choices[0].message.content
        except Exception as e:
            retries += 1
            print(f"Error occurred, retrying ({retries}/{max_retries}): {e}")
            time.sleep(retry_delay)
    
    return "Max retries exceeded, unable to get a response from the API."


def sys_question(sys_prompt, prompt, max_retries=5, retry_delay=5):
    """
    添加系统提示词的API调用, print(ask_gpt_question("无论我输入什么，你都输出你好","1+1=?")
    """
    retries = 0
    while retries < max_retries:
        try:
            response = openai.chat.completions.create(
                model='o4-mini-high-all',  
                messages=[
                    {"role":"system","content":sys_prompt},
                    {"role": "user", "content": prompt}
                ],
            )
            return response.choices[0].message.content
        except Exception as e:
            retries += 1
            print(f"Error occurred, retrying ({retries}/{max_retries}): {e}")
            time.sleep(retry_delay)
    
    return "Max retries exceeded, unable to get a response from the API."


def patch_sys_question(sys_prompt, prompt, before_file, max_retries=1, retry_delay=5):
    """
    基于系统提示词API改进的patch生成, 用于生成可执行patch, 示例见draft/test_API_patch.py
    """
    retries = 0
    while retries < max_retries:
        try:
            response = openai.chat.completions.create(
                model='o4-mini-high-all',
                messages=[
                    {"role": "system", "content": sys_prompt},
                    {"role": "user", "content": prompt}
                ],
            )
            api_response = response.choices[0].message.content
            patch_content = re.search(r'```patch(.*?)```', api_response, re.DOTALL)
            retries += 1
            if patch_content:
                patch_content = patch_content.group(1).strip()
                patch_content = patch_content + '\n'
                with open(f'{base_dir}/draft/test.patch', 'w') as f:
                    f.write(patch_content)

                try:
                    subprocess.run(['patch', before_file, '-o', f'{base_dir}/draft/after.c'], input=patch_content.encode(), check=True)
                    return patch_content
                except subprocess.CalledProcessError as e:
                    print(f"Error applying patch: {e}")
            else:
                print("No patch content found in the API response.")

        except Exception as e:
            retries += 1
            print(f"Error occurred, retrying ({retries}/{max_retries}): {e}")
            print('*'*10)
            time.sleep(retry_delay)

    return "Max retries exceeded, unable to get a response from the API."


async def ask_gpt_question_async(prompt, max_retries=5, retry_delay=5):
    retries = 0
    while retries < max_retries:
        try:
            response = await client.chat.completions.create(
                model='gpt-4o-mini',
                messages=[{"role": "user", "content": prompt}]
            )
            return response.choices[0].message.content
        except Exception as e:
            retries += 1
            print(f"Error occurred, retrying ({retries}/{max_retries}): {e}")
            await asyncio.sleep(retry_delay)
    return "Max retries exceeded, unable to get a response from the API."

async def batch_ask_questions(questions):
    tasks = [ask_gpt_question_async(q) for q in questions]
    return await asyncio.gather(*tasks)

