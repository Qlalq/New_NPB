import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from API import ask_gpt_question, patch_sys_question

base_dir = os.environ.get("BASE_DIR")
before_address = f"{base_dir}/draft/add.c"
sys_address = f"{base_dir}/prompt/0.txt"
with open(sys_address, 'r') as f:
    sys_prompt = f.read()

prompt="""
static void add(void) {
  int i, j, k, m;
  for (i = 1; i < grid_points[0]-1; i++) {
    for (j = 1; j < grid_points[1]-1; j++) {
      for (k = 1; k < grid_points[2]-1; k++) {
	for (m = 0; m < 5; m++) {
	  u[i][j][k][m] = u[i][j][k][m] + rhs[i][j][k][m];
	}
      }
    }
  }
}
"""


print(patch_sys_question(sys_prompt, prompt, before_address))