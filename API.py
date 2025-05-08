import openai
import time
import re
import os
import subprocess
base_dir = os.environ.get("BASE_DIR")

openai.base_url = 'https://api.siliconflow.cn/v1/'
openai.api_key = ''

def ask_gpt_question(prompt, max_retries=5, retry_delay=5):
    """
    最简单的API调用, print(ask_gpt_question("hello")
    """
    retries = 0
    while retries < max_retries:
        try:
            response = openai.chat.completions.create(
                model='gemini-2.5-flash-preview-04-17',  
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
                model='gemini-2.5-flash-preview-04-17',  
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



def patch_sys_question(sys_prompt, prompt, before_file, max_retries=3, retry_delay=5):
    """
    基于系统提示词API改进的patch生成, 用于生成可执行patch, 示例见draft/test_API_patch.py
    """
    retries = 0
    while retries < max_retries:
        try:
            response = openai.chat.completions.create(
                model='gemini-2.5-flash-preview-04-17',
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
            time.sleep(retry_delay)

    return "Max retries exceeded, unable to get a response from the API."
