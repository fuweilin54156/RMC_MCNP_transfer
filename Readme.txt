MCNP to RMC转换工具 （V2.0）
开发者：申鹏飞，苟元昊
邮箱：2043965149@qq.com

此脚本可以转换MCNP输入卡中的几何、材料、重要性、TR卡、源描述卡等部分，并生成对应的RMC输入卡。

使用说明：
1. MCNPtoRMC.py 中提供了转换函数，可参考 runner.py 文件调用并处理。可以通过运行 runner.py 脚本或 M2Rtransfer 程序的方式开始转换，输入文件名，自动调用转换函数处理。支持‘？’和‘*’的通配输入（仅一次）。
其中，‘？’匹配单个字符，‘*’匹配任意字符。
在命令行的转换界面，会以 Warning: 和 Error: 的形式提示未转换部分，请务必仔细核对！

2. 针对每个输入的MCNP文件，会输出三个文件：
（1）[输入文件名].rmc.python
已完成解析的RMC的Python包可读取的文件，如果不含有RMC宏体等部分，也可以由RMC二进制程序直接读取
（2）[输入文件名].rmc.binary
已完成解析的RMC的二进制程序可读取的文件
（3）[输入文件名].mcnp
已完成解析的MCNP的模型文件，会以 Warning: 和 Error: 的形式提示未转换部分，请务必仔细核对！

注意：
1. MCNP的输入卡表示Cell和Surf卡每一行的数字编号前不能有空格（否则会转换失败）
2. 转换后需要重点检查RMC的lattice参数、tr参数

参考资料：
RMC用户手册在线版：https://rmc-doc.reallab.org.cn/en/latest/usersguide/index.html
