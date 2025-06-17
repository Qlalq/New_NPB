from project import init_NPB, replace_NPB, run_NPB

# 设置参数
benchmark = "CG"
function_name = "conj_grad"  # 举例
folder_name = "function_refinement"
class_type = "S"  # 可选 S, W, A, B, C 等

# 初始化原始源文件
init_NPB(benchmark)

# 替换函数为你生成的新版本
# replace_NPB(benchmark, function_name, folder_name)

# 编译并运行 NPB，返回正确性与耗时
success, runtime = run_NPB(benchmark, class_type)

if success:
    print(f"[✓] 验证成功，运行时间：{runtime:.3f} 秒")
else:
    print("[✗] 验证失败或运行异常")