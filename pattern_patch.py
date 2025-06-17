import os
import re
import glob
import subprocess
import argparse
from API import patch_sys_question
from project import code_line

def process_folder(base_dir, folder_name, prompt_file):
    npb_base_folder = f'{base_dir}/NPB3.0-omp-C'
    source_folder = f'{npb_base_folder}/{folder_name}/function_#_omp'
    baseline_folder = f'{npb_base_folder}/{folder_name}/function_pattern_baseline'
    patch_folder = f'{baseline_folder}/patch'
    output_folder = f'{baseline_folder}/C_patch'

    # 确保目标文件夹存在
    os.makedirs(baseline_folder, exist_ok=True)
    os.makedirs(patch_folder, exist_ok=True)
    os.makedirs(output_folder, exist_ok=True)

    # 查找源文件夹中的所有 C 文件
    c_files = glob.glob(os.path.join(source_folder, "*.c"))

    for c_file_path in c_files:
        # 获取文件名
        c_file_name = os.path.basename(c_file_path)

        # 读取 C 文件内容并获取行号
        codeline_input = code_line(c_file_path)
        # 生成 patch
        print(f"Processing {c_file_name} in {folder_name}...")
        patch_content = patch_sys_question(
            prompt_file,
            codeline_input,
            c_file_path
        )

        # 将 patch 保存到文件
        patch_file_path = os.path.join(patch_folder, f"{c_file_name}.patch")
        with open(patch_file_path, 'w') as f:
            f.write(patch_content)
        print(f"Saved patch file to {patch_file_path}")

        # 应用 patch 并保存修改后的文件
        output_file_path = os.path.join(output_folder, c_file_name)
        try:
            subprocess.run(['patch', c_file_path, '-o', output_file_path], input=patch_content.encode(), check=True)
            print(f"Saved refined C file to {output_file_path}")
        except subprocess.CalledProcessError as e:
            print(f"Error applying patch: {e}")

    print(f"Finished processing folder: {folder_name}")


def main():
    parser = argparse.ArgumentParser(description="pattern baseline generation")
    parser.add_argument(
        '--folder', 
        nargs='+', 
        default=['BT'],
        help='List of folder names to process (e.g., BT CG SP)'
    )
    args = parser.parse_args()

    base_dir = os.environ.get("BASE_DIR")
    prompt_file = f"{base_dir}/prompt/0_patch.txt"
    with open(prompt_file, 'r') as f:
        prompt_template = f.read()

    # 遍历所有指定的文件夹名称
    for folder_name in args.folder:
        process_folder(base_dir, folder_name, prompt_template)

    print("All folders processed!")


if __name__ == "__main__":
    main()