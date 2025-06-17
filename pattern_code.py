import os
import glob
import re
import argparse
from API import sys_question
from project import code_line

def extract_code_from_markdown(generated_text):
    """Extract only the code content from markdown code blocks."""
    pattern = r"```c\s*(.*?)\s*```"
    matches = re.findall(pattern, generated_text, re.DOTALL)
    
    if matches:
        return matches[0]  # Return the first code block found
    else:
        return generated_text

def process_folder(base_dir, folder_name):
    npb_base_folder = f'{base_dir}/NPB3.0-omp-C'
    source_folder = f'{npb_base_folder}/{folder_name}/function_#_omp'
    baseline_folder = f'{npb_base_folder}/{folder_name}/function_pattern_baseline'
    output_folder = f'{baseline_folder}/C_code'  # Updated output folder

    # Ensure target folders exist
    os.makedirs(baseline_folder, exist_ok=True)
    os.makedirs(output_folder, exist_ok=True)

    # Find all C files in the source folder
    c_files = glob.glob(os.path.join(source_folder, "*.c"))

    for c_file_path in c_files:
        # Get the file name
        c_file_name = os.path.basename(c_file_path)

        # Read C file content and get code lines
        codeline_input = code_line(c_file_path)
        print(f"Processing {c_file_name} in {folder_name}...")

        # Generate C code directly using the sys_question API
        prompt_file = f"{base_dir}/prompt/0_code.txt"
        with open(prompt_file, 'r') as f:
            sys_prompt = f.read()

        # Call the API to generate C code
        generated_response = sys_question(sys_prompt, codeline_input)
        
        # Extract only the C code from the response
        filtered_code = extract_code_from_markdown(generated_response)

        # Save the generated C code to the output folder
        output_file_path = os.path.join(output_folder, c_file_name)
        with open(output_file_path, 'w') as f:
            f.write(filtered_code)
        print(f"Saved generated C file to {output_file_path}")

    print(f"Finished processing folder: {folder_name}")

def main():
    parser = argparse.ArgumentParser(description="Pattern baseline generation")
    parser.add_argument(
        '--folder', 
        nargs='+', 
        default=['BT'],
        help='List of folder names to process (e.g., BT CG SP)'
    )
    args = parser.parse_args()

    base_dir = os.environ.get("BASE_DIR")
    
    # Iterate over all specified folder names
    for folder_name in args.folder:
        process_folder(base_dir, folder_name)

    print("All folders processed!")

if __name__ == "__main__":
    main()