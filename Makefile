#----------------------------------------------------#
#--Author: Harrys Kon (Charalambos Konstantinou)-----#
#--W: https://harrys.fyi/----------------------------#
#--E: konharrys@gmail.com----------------------------#
#----------------------------------------------------#
#CPP = g++
DEBUG = -g
CPP_ANSI_OPTS = -std=c++14
RM=rm -f
EIGEN = `pkg-config --libs --cflags eigen3`
EIGEN1 = -I /usr/local/include/eigen3
CPPFLAGS = -Wall $(DEBUG) $(CPP_ANSI_OPTS)
SRC_DIR = src

ifeq "$(shell expr `g++ -dumpversion | cut -f1 -d.` \> 4.9)" "1"
	CPP = g++
else ifeq "$(shell expr `g++-4.9 -dumpversion | cut -c 1-3` \= 4.9)" "1"
	CPP = g++-4.9
else
	$(error Please install gcc/g++ 4.9 or greater)
endif


SOURCES= $(wildcard $(SRC_DIR)/*.cpp)
EXECUTABLE=nrpowerflow
    
$(EXECUTABLE): 
	$(CPP) $(CPPFLAGS) $(LDFLAGS) $(SOURCES) -o $@ $(EIGEN)

.PHONY : $(EXECUTABLE) clean

clean :
	$(RM) $(EXECUTABLE)
