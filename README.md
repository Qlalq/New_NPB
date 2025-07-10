# NPB并行化测试项目

## 项目简介

这是一个项目级的并行框架，提供了在NPB (NAS Parallel Benchmarks) 上优化的初步尝试，旨在探索通过自动化手段（如代码重构和针对性优化）来提升C/C++项目的并行性能。

项目主要包含以下两个核心功能：

1.  **代码重构**
2.  **针对性优化**

这两个功能分别由项目根目录下的 `refine.py` 和 `pattern.py` 脚本实现。

## 快速开始

### 1. 环境变量设置

在运行项目脚本之前，需要设置 `BASE_DIR`和`API_KEY` 环境变量，指向项目的根目录（即 `New_NPB_frame` 所在的路径）。

```shell
export BASE_DIR="/path/to/your/New_NPB_frame"
# 例如：
# export BASE_DIR="workspace/New_NPB_frame"
```

### 2. 配置 API Key

请在 `API.py` 文件中填入项目所需的API Key，以便调用相关的服务（如大型语言模型等）。

```python
# API.py 文件示例 (请替换 your_api_key_here)
# ... 其他导入或代码 ...
openai.api_key = "your_api_key_here"
# ... 其他API调用代码 ...
```


## 使用方法

### 快速开始

*   **进行代码重构：**

    ```bash
    python refine.py
    ```

*   **进行针对性优化：**

    ```bash
    python pattern.py
    ```

## 文件结构说明

项目根目录包含以下主要文件和文件夹：

-   `API.py`: 负责处理与API相关的各类调用，需要在此文件中配置API Key。
-   `project.py`: 封装了项目级操作的相关逻辑。
### `draft/`
存储在项目运行过程中产生的临时文件。

### `NPB3.0-omp-C/`

该目录下存放NPB基准测试的原始代码、处理后的代码和相关数据。它分为7大类基准测试：

-   BT
-   CG
-   EP
-   FT
-   LU
-   MG
-   SP

每个类别（例如 `BT/`）的文件和文件夹结构类似，下面以 `BT/` 目录为例进行讲解：

-   `bt_ori.c`: 人工专家优化后的原始C语言源代码备份。
-   `bt_#_omp.c`: 基于 `bt_ori.c` 去除注释和OpenMP原语后的备份。
-   `bt.c`: 实际用于编译和执行测试的源代码文件。
-   `extra_fun.py`: 用于从源代码文件中提取独立函数的功能脚本。
-   `function_#_omp/`: 存储从 `bt_#_omp.c` 中提取的独立函数文件，作为重构雨针对优化的数据。
-   `function_refinement/`: 存储代码重构功能（`refine.py`）生成的重构后的函数文件。现有的文件是基于Gemini2.5Flash生成的，并非完全正确，具体见`NPB3.0-omp-C/All/refine.jsonl`。
-   `function_pattern_baseline/`: 存储针对性优化功能（`pattern.py`）生成的优化后的函数文件，其下有两个文件夹，`C/`存储优化后的c文件，`patch/`存储LLM生成的patch。现有的文件是基于Gemini2.5Flash生成的，并非完全正确，空白的C文件表示生成失败。

### `prompt/`

存放项目系统使用的提示词（Prompt）文件，用于指导AI模型进行代码处理。

-   `0.txt`: 通用优化提示词。
-   `1.txt` 到 `5.txt`: 针对性优化的不同提示词，可能对应不同的优化策略或模式。
-   `refine.txt` : 代码重构功能所使用的提示词文件。

