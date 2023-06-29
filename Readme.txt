MCNP to RMC转换工具 （V2.0）
开发者：申鹏飞，苟元昊
邮箱：2043965149@qq.com

使用说明：
1. MCNPtoRMC.py 中提供了转换函数，可参考 runner.py 文件调用并处理。运行 runner.py 可以输入文件名，自动调用转换函数处理。支持‘？’和‘*’的通配输入（仅一次）。
其中，‘？’匹配单个字符，‘*’匹配任意字符
2. 针对每个输入的MCNP文件，会输出三个文件：
（1）[输入文件名].rmc.python
已完成解析的RMC的Python包可读取的文件，如果不含有RMC宏体等部分，也可以由RMC二进制程序直接读取
（2）[输入文件名].rmc.binary
已完成解析的RMC的二进制程序可读取的文件
（3）[输入文件名].mcnp
已完成解析的MCNP的模型文件，会以 warning： 的形式提示未转换部分，请务必仔细核对

注意：
1. MCNP的输入卡表示Cell和Surf的每一行 数字前 不能有空格（否则会转换失败）
2. 转换后需要重点检查RMC的lattice参数

参考资料：
RMC用户手册在线版：https://rmc-doc.reallab.org.cn/en/latest/usersguide/index.html
