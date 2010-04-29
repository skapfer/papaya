import sys

class ABORT:
    pass

def die (msg):
    print >>sys.stderr, msg
    raise ABORT ()


# ***************************************************************************
# break line into strings
#
# SPECIAL VERSION 
#
# this function breaks a string 'line' as read by
# line=file.readline()
# into separate strings where it recognises blanc spaces, '(', ')' and  ':' as
# separators between these strings.
# "asdfsd 234 dfsdf   234" -> {"asdfsd","234","dfsdf","234"}
#
# A line read from a file probably ends with a carriage return,
# this will not be returned
#
def break_line_into_string(line):
    list=[]
    length=len(line)
    i=0
    if line[len(line)-1] == '\012':
        line = line[:-1]
    while i < len(line):
        while i<len(line):
            break_now = 1
            if line[i]==' ' or line[i] == ':' or line[i] == '\t' or \
               line[i]=='(' or line[i] == ')' or line[i] == ',' or\
               line[i]=='[' or line[i] == ']' or line[i] == '=':
                break_now = 0
                i=i+1
            if break_now == 1:
                break
        word=""
        while i<len(line):
            if line[i] == ' ' or line[i] == ':'  or line[i] == '\t' or\
               line[i] == '(' or line[i] == ')' or line[i] == ',' or\
               line[i]=='[' or line[i] == ']' or line[i] == '=':
                break
            word=word+line[i]
            i=i+1
        if word != "":
            list.append(word)
    return list
