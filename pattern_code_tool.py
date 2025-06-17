import os
import glob
import re
import argparse
import subprocess
import shutil
from API import *
from project import *

def extract_code_from_markdown(generated_text):
    """Extract only the code content from markdown code blocks."""
    pattern = r"```c\s*(.*?)\s*```"
    matches = re.findall(pattern, generated_text, re.DOTALL)
    
    if matches:
        return matches[0]  # Return the first code block found
    else:
        return generated_text

def extract_omp_primitives(response, function_name):
    """Extract OMP primitives for a specific function from response."""
    pattern = f"```{function_name}\n(.*?)```"
    match = re.search(pattern, response, re.DOTALL)
    
    if match:
        return match.group(1).strip()
    return ""

def backup_file(file_path):
    """Create a backup of a file."""
    backup_path = f"{file_path}_tmp"
    shutil.copy2(file_path, backup_path)
    return backup_path

def restore_file(backup_path, dest_path):
    """Restore a file from its backup."""
    shutil.copy2(backup_path, dest_path)
    return True

def process_folder(base_dir, folder_name, CLASS="S", retry=3):
    npb_base_folder = f'{base_dir}/NPB3.0-omp-C'
    source_folder = f'{npb_base_folder}/{folder_name}/function_#_omp'
    baseline_folder = f'{npb_base_folder}/{folder_name}/function_pattern_baseline'
    output_folder = f'{baseline_folder}/C_code'  # Updated output folder
    init_NPB(folder_name)
    # Ensure target folders exist
    os.makedirs(baseline_folder, exist_ok=True)
    os.makedirs(output_folder, exist_ok=True)

    bench_lower = folder_name.lower()
    bench_file = f"{npb_base_folder}/{folder_name}/{bench_lower}.c"
    
    # Initialize a dictionary to store function OMP primitives
    function_omp = {}
    
    success, run_time = run_NPB(folder_name, CLASS)
    init_time = run_time
    print(f'init_time:{init_time}')
    
    # 首先分析主程序文件以确定函数调用关系和BFS顺序
    functions, dependencies = analyze_c_function_calls(bench_file)
    
    # Get the BFS order for optimization
    bfs_order = get_bfs_order(functions, dependencies)
    
    print(f"Optimization order for functions in {folder_name}: {bfs_order}")
    
    # Process functions in the optimized order
    for function_name in bfs_order:
        # 在source_folder中查找对应函数的源文件
        function_file_path = os.path.join(source_folder, f"{function_name}.c")
        
        if not os.path.exists(function_file_path):
            print(f"Warning: File for function {function_name} not found in {source_folder}, skipping...")
            continue
            
        print(f"Processing {function_name}.c in {folder_name}...")
        
        # Read C file content and get code lines
        with open(function_file_path, 'r') as f:
            codeline_input = f.read()
        
        # Retry loop for code generation
        success_optimization = False
        for attempt in range(retry):
            print(f"Attempt {attempt+1}/{retry} for optimizing {function_name}")
            
            # Generate C code using the sys_question API
            prompt_file = f"{base_dir}/prompt/0_code.txt"
            with open(prompt_file, 'r') as f:
                sys_prompt = f.read()
            
            # Add OMP primitives information to the prompt if available
            if function_name in function_omp and function_omp[function_name]:
                omp_info = f"{function_omp[function_name]}"
                sys_prompt += omp_info
            
            # Call the API to generate C code
            generated_response = sys_question(sys_prompt, codeline_input)
            
            # Extract only the C code from the response
            filtered_code = extract_code_from_markdown(generated_response)
            print(f'\nGeneration Code:\n{filtered_code} ')

            # Save the generated C code to the output folder
            output_file_path = os.path.join(output_folder, f"{function_name}.c")
            with open(output_file_path, 'w') as f:
                f.write(filtered_code)
            print(f"Saved generated C file to {output_file_path}")
            
            # Create a backup of the original file before replacement
            backup_path = backup_file(bench_file)
            print(f"Created backup of {bench_file} to {backup_path}")
            
            # Replace in NPB and check performance
            replace_NPB(folder_name, function_name, f"function_pattern_baseline/C_code")
            success, run_time = run_NPB(folder_name, CLASS)
            print('*'*10)
            print(f'init_time:{init_time}')
            print(f'best_time:{best_time}')
            print(f'{success},{run_time}')
            print('*'*10)
            
            if success and run_time < best_time:
                best_time = run_time
                print(f"New best time for {folder_name}: {best_time} seconds with optimization of {function_name}")
                
                # After successful replacement, analyze OMP context for called functions
                if function_name in dependencies and dependencies[function_name]:
                    called_functions = dependencies[function_name]
                    
                    # Prepare prompt for finding OMP primitives
                    with open(f"{base_dir}/prompt/find_omp.txt", 'r') as f:
                        sys_prompt_context = f.read()
                    
                    # Add the list of called functions
                    extra_content = ", ".join(called_functions)
                    sys_prompt_context += extra_content
                    
                    # Get the updated code
                    with open(output_file_path, 'r') as f:
                        updated_code = f.read()
                    
                    # Call API to find OMP primitives affecting called functions
                    omp_response = sys_question(sys_prompt_context, updated_code)
                    
                    # Extract and store OMP primitives for each called function
                    for called_func in called_functions:
                        omp_primitives = extract_omp_primitives(omp_response, called_func)
                        if omp_primitives:
                            function_omp[called_func] = omp_primitives
                            print(f"Found OMP primitives for {called_func}: {omp_primitives}")
                
                # Delete backup if successful
                os.remove(backup_path)
                success_optimization = True
                break  # Exit retry loop on success
            else:
                # If not successful or not better, restore from backup
                print(f"Attempt {attempt+1} failed. No improvement or execution failed.")
                restore_file(backup_path, bench_file)
                print(f"Restored {bench_file} from backup.")
                # Delete backup after restore
                os.remove(backup_path)
                
                # Continue to next retry attempt if we haven't exhausted retries
                if attempt < retry - 1:
                    print(f"Retrying optimization for {function_name}...")
                else:
                    print(f"All {retry} attempts for {function_name} failed to improve performance.")
        
        if not success_optimization:
            print(f"Failed to optimize {function_name} after {retry} attempts. Moving to next function.")

    print(f"Finished processing folder: {folder_name}")
    print(f"Best execution time: {best_time} seconds")
    print(f"OMP primitives collected: {function_omp}")


def main():
    parser = argparse.ArgumentParser(description="Pattern baseline generation with optimization ordering")
    parser.add_argument(
        '--folder', 
        nargs='+', 
        default=['CG'],
        help='List of folder names to process (e.g., BT CG SP)'
    )
    parser.add_argument(
        '--class', 
        default='W',
        help='Problem class size (S, W, A, B, C)'
    )
    args = parser.parse_args()

    base_dir = os.environ.get("BASE_DIR")
    
    # Iterate over all specified folder names
    for folder_name in args.folder:
        process_folder(base_dir, folder_name, getattr(args, 'class'))

    print("All folders processed!")

if __name__ == "__main__":
    main()