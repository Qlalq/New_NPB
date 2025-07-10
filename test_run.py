from project import *

# 设置参数
# benchmarks = ["BT", "CG", "EP", "FT", "LU", "MG", "SP"]
benchmarks = ["CG", "MG"]
class_type = "W"  # 可选 S, W, A, B, C 等
num_runs = 1  # 运行次数

# 存储所有结果
results = {}

for benchmark in benchmarks:
    print(f"\n正在测试 {benchmark}...")
    
    # 初始化原始源文件
    init_NPB(benchmark)
    
    # 记录所有运行时间
    runtimes = []
    all_success = True
    
    for run_idx in range(num_runs):
        print(f"  第 {run_idx + 1}/{num_runs} 次运行...")
        success, runtime = run_NPB(benchmark, class_type,times=3)
        
        if not success:
            print(f"[✗] {benchmark} 验证失败或运行异常，跳过此benchmark")
            all_success = False
            break
        else:
            runtimes.append(runtime)
            print(f"    运行时间：{runtime:.3f} 秒")
    
    # 存储结果
    if all_success:
        avg_runtime = sum(runtimes) / len(runtimes)
        results[benchmark] = {
            'success': True,
            'avg_runtime': avg_runtime,
            'runtimes': runtimes
        }
    else:
        results[benchmark] = {
            'success': False,
            'avg_runtime': None,
            'runtimes': []
        }

print("\n" + "="*50)
print("测试结果汇总：")
print("="*50)

for benchmark in benchmarks:
    if results[benchmark]['success']:
        avg_time = results[benchmark]['avg_runtime']
        detail_times = results[benchmark]['runtimes']
        print(f"[✓] {benchmark} 验证成功，平均运行时间：{avg_time:.3f} 秒")
        print(f"    详细时间：{[f'{t:.3f}' for t in detail_times]} 秒")
    else:
        print(f"[✗] {benchmark} 验证失败或运行异常")

print("="*50)
print("所有benchmark测试完成")