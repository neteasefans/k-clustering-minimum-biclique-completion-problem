#以下是指定编译器路径(Here are the specified compiler paths)
CC = /usr/bin/g++
#以下是指定编译需要的头文件(Here are the header files required for compilation.)
CFLAGS = -g -Wall -o3
#以下是源文件(Here are the .cpp files)
SRCS = ./code/*.cpp
#以下是指定目标文件 所有的.c文件变成.o文件
OBJS = $(SRCS:.cpp=.o)
#以下是生成可执行文件(The executable file)
EXECUTABLE = main_exe


#make all 执行生成可执行文件
#1编译器 2编译选项 3输出 4生成的可执行文件 5需要的源文件 6需要的库文件
all:
	$(CC) $(CFLAGS) -o $(EXECUTABLE) $(SRCS) 
clean:
	rm -f *.o 
