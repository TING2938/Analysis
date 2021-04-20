# 计算ECPM矩阵

## 1. 安装

```bash
tar -zxvf calcMatrix.tar.gz
cd calcMatrix
chmod +x configure INSTALL
./INSTALL EcpmMatrix/calcMatrix.cpp -DUSE_THREAD=ON
```

会在`./bin`目录下生成可执行文件`calcMatrix`.

## 2.使用方法

- 输出帮助信息

  ```bash
  $ ./calcMatrix -h 
  Command line option:
  calcMatrix -h
  
  Sub Program Message:
         Function     Description
  --------------------------------------------------------
           matrix     Calculate Matrix for ECPM
          potFile     Get "getPot_parameters.dat" for different voltage
  --------------------------------------------------------
  ```

  ```bash
  $ ./calcMatrix matrix -h
  calculate Matrix for ECPM
  
  Command line option:
  calcMatrix matrix -h
  
           Option                    Value     Description
  --------------------------------------------------------
  (Required)
               -f                 wall.gro     gro file name
               -n                        0     number of total electrode atoms
              -kx                        1     kxmax
              -ky                        1     kymax
              -kz                        1     kzmax
             -3dc                        1     3d(0) or 3dc(1)
  (Optional)
          -thread                        8     number of thread
          -cutoff                      1.2     cut off(nm)
             -eta                    1.979     eta
           -ewald                 0.260284     g_ewald
  --------------------------------------------------------
  ```

  ```bash
  $ ./calcMatrix potFile -h
  get "getPotFile" for ECPM
  
  Command line option:
  calcMatrix potFile -h
  
           Option                    Value     Description
  --------------------------------------------------------
  (Required)
               -f                 wall.gro     gro file name
               -n                        0     number of total electrode atoms
               -N                        0     number of total atoms of system
               -v                [0, 1, 2]     voltage
             -low                      0.0     bulk low position(z)(nm)
              -up                     10.0     bulk up position(z)(nm)
  (Optional)
          -cutoff                      1.2     cut off(nm)
  --------------------------------------------------------
  ```
  
- 示例

  - 计算矩阵：

    ```bash
    ./calcMatrix matrix -f wall.gro -n 2178 -kx 20 -ky 20 -kz 30 -3dc 0
    ```

  - 计算`getPot_parameters.dat`文件

    ```bash
    ./calcMatrix potFile -f wall.gro -n 2178 -N 24700 -v 0 0.5 1 1.5 2 -low 1.5 -up 3.5
    ```

ting ye

yeting2938@hust.edu.cn

07/30 2020