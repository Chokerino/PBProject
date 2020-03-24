#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<string.h>

void executeScript(char* GEOID)
{
    
    pid_t pid = fork();
    if(pid>0){
        char printPID[1024] = "echo Process PID: ";
        char processID[24];
        sprintf(processID, "%d", pid);
        strcat(printPID,processID);
        strcat(printPID," >> log.txt");
        system(printPID);
    }
    if(pid<=0)
    {
        char command[1024] = "python3 nvbi.py ";
        strcat(command,GEOID);
        strcat(command," >> log.txt");
        system(command);
    }
    return;
}

int main(int argn, char** argc)
{
    if(argn<=1)
    {
        printf("Geoid not present\n");
        return -1;
    }
    executeScript(argc[1]);
    return 0;
}