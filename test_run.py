from project import *
import numpy as np
import matplotlib.pyplot as plt

# 设置参数
# benchmarks = ["BT", "CG", "EP", "FT", "LU", "MG", "SP"]
benchmarks = ["EP"]
class_type = "S"  # 可选 S, W, A, B, C 等
num_runs = 150  # 运行次数

# 存储所有结果
results = {}

for benchmark in benchmarks:
    print(f"\n正在测试 {benchmark}...")
    
    # # 初始化原始源文件
    # init_NPB(benchmark, "autoPar")
    
    # 记录所有运行时间
    runtimes = []
    cumulative_means = []  # 累计平均值
    cumulative_vars = []   # 累计方差
    all_success = True
    
    for run_idx in range(num_runs):
        print(f"  第 {run_idx + 1}/{num_runs} 次运行...")
        success, runtime = run_NPB(benchmark, class_type)
        
        if not success:
            print(f"[✗] {benchmark} 验证失败或运行异常，跳过此benchmark")
            all_success = False
            break
        else:
            runtimes.append(runtime)
            print(f"    运行时间：{runtime:.3f} 秒")
            
            # 计算累计平均值和方差
            current_mean = np.mean(runtimes)
            current_var = np.var(runtimes, ddof=1) if len(runtimes) > 1 else 0
            
            cumulative_means.append(current_mean)
            cumulative_vars.append(current_var)
            
            print(f"    累计平均：{current_mean:.3f} 秒，方差：{current_var:.6f}")
    
    # 存储结果
    if all_success:
        final_avg = np.mean(runtimes)
        final_var = np.var(runtimes, ddof=1)
        final_std = np.std(runtimes, ddof=1)
        
        results[benchmark] = {
            'success': True,
            'avg_runtime': final_avg,
            'variance': final_var,
            'std_dev': final_std,
            'runtimes': runtimes,
            'cumulative_means': cumulative_means,
            'cumulative_vars': cumulative_vars
        }
    else:
        results[benchmark] = {
            'success': False,
            'avg_runtime': None,
            'variance': None,
            'std_dev': None,
            'runtimes': [],
            'cumulative_means': [],
            'cumulative_vars': []
        }

print("\n" + "="*50)
print("测试结果汇总：")
print("="*50)

for benchmark in benchmarks:
    if results[benchmark]['success']:
        avg_time = results[benchmark]['avg_runtime']
        variance = results[benchmark]['variance']
        std_dev = results[benchmark]['std_dev']
        detail_times = results[benchmark]['runtimes']
        
        print(f"[✓] {benchmark} 验证成功")
        print(f"    平均运行时间：{avg_time:.3f} 秒")
        print(f"    方差：{variance:.6f}")
        print(f"    标准差：{std_dev:.6f} 秒")
        print(f"    变异系数：{(std_dev/avg_time)*100:.2f}%")
        print(f"    最小时间：{min(detail_times):.3f} 秒")
        print(f"    最大时间：{max(detail_times):.3f} 秒")
    else:
        print(f"[✗] {benchmark} 验证失败或运行异常")

# 生成图表
plt.rcParams['axes.unicode_minus'] = False

for benchmark in benchmarks:
    if results[benchmark]['success']:
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 10))
        
        iterations = range(1, len(results[benchmark]['cumulative_means']) + 1)
        cumulative_means = results[benchmark]['cumulative_means']
        cumulative_vars = results[benchmark]['cumulative_vars']
        runtimes = results[benchmark]['runtimes']
        
        # 子图1：累计平均值变化
        ax1.plot(iterations, cumulative_means, 'b-', linewidth=2, label='Cumulative average value')
        ax1.axhline(y=results[benchmark]['avg_runtime'], color='r', linestyle='--', 
                   label=f'Final average value: {results[benchmark]["avg_runtime"]:.3f}s')
        ax1.set_xlabel('Number of iterations')
        ax1.set_ylabel('Average running time (s)')
        ax1.set_title(f'{benchmark} - The convergence process of cumulative average values')
        ax1.grid(True, alpha=0.3)
        ax1.legend()
        
        # 子图2：累计方差变化
        ax2.plot(iterations, cumulative_vars, 'g-', linewidth=2, label='Cumulative variance')
        ax2.axhline(y=results[benchmark]['variance'], color='r', linestyle='--', 
                   label=f'Final variance: {results[benchmark]["variance"]:.6f}')
        ax2.set_xlabel('Number of iterations')
        ax2.set_ylabel('Variance')
        ax2.set_title(f'{benchmark} - Cumulative variance change')
        ax2.grid(True, alpha=0.3)
        ax2.legend()
        
        # 子图3：每次运行时间散点图
        ax3.scatter(iterations, runtimes, alpha=0.6, s=20, c='orange', label='Single run time')
        ax3.axhline(y=results[benchmark]['avg_runtime'], color='r', linestyle='-', 
                   alpha=0.8, label=f'Average value: {results[benchmark]["avg_runtime"]:.3f}s')
        # 添加±1σ和±2σ区间
        mean_val = results[benchmark]['avg_runtime']
        std_val = results[benchmark]['std_dev']
        ax3.axhline(y=mean_val + std_val, color='orange', linestyle=':', alpha=0.7, label='±1σ')
        ax3.axhline(y=mean_val - std_val, color='orange', linestyle=':', alpha=0.7)
        ax3.axhline(y=mean_val + 2*std_val, color='red', linestyle=':', alpha=0.7, label='±2σ')
        ax3.axhline(y=mean_val - 2*std_val, color='red', linestyle=':', alpha=0.7)
        
        ax3.set_xlabel('Number of iterations')
        ax3.set_ylabel('Average running time (s)')
        ax3.set_title(f'{benchmark} - Single-run time distribution')
        ax3.grid(True, alpha=0.3)
        ax3.legend()
        
        plt.tight_layout()
        plt.savefig(f'{benchmark}_runtime_analysis.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        # 打印收敛性分析
        print(f"\n{benchmark} 收敛性分析：")
        # 计算最后10%的数据的稳定性
        last_10_percent = int(num_runs * 0.1)
        if last_10_percent > 0:
            recent_means = cumulative_means[-last_10_percent:]
            recent_mean_var = np.var(recent_means)
            print(f"  最后{last_10_percent}次迭代的平均值方差：{recent_mean_var:.8f}")
            print(f"  平均值变化趋势：{'收敛' if recent_mean_var < 0.001 else '仍在波动'}")

print("="*50)
print("所有benchmark测试完成")