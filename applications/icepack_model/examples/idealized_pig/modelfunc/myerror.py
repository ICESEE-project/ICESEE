import sys
def myerror(message,myLogger=None):
    """ print error and exit """
    print('\n\t\033[1;31m *** ',message,' *** \033[0m\n') 
    if myLogger != None :
        myLogger.logError(message)
    sys.exit()




