# PhysicsEducation

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://wwangnju.github.io/PhysicsEducation.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://wwangnju.github.io/PhysicsEducation.jl/dev/)
[![Build Status](https://github.com/wwangnju/PhysicsEducation.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/wwangnju/PhysicsEducation.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/wwangnju/PhysicsEducation.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/wwangnju/PhysicsEducation.jl)

## 安装
In Julia, please type `]` in the REPL to use the package mode, then type this command:

```julia
pkg> add https://github.com/wwangnju/PhysicsEducation.jl.git
```

## 例子
在examples目录中，实现了光栅夫琅和费衍射、圆孔夫琅和费衍射、简谐振动、受迫振动、牛顿环、双缝干涉、驻波以及正多边形夫琅和费衍射的可交互式动态模拟。
![standing](./examples/驻波.png)
每次改变参数后，点击“更新”按钮。“开始”按钮点击后出现动画演示，再次点击会暂停。“初始化”按钮将轨迹初始化到更新参数后的位置。“前进”按钮点击可以实现暂停后动画前进一步，“步长”表示每次前进的轨迹长度。
![niu](./examples/牛顿环.png)
