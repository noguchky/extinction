# coding: utf-8

import sys
import cmd
import os
import struct
import sitcprbcp
import fileinput

DEFAULT_IP = '192.168.10.16'
DEFAULT_PORT = '4660'

rbcpId = 0;

class rbcpCmd(cmd.Cmd):
    def __init__(self):
        cmd.Cmd.__init__(self)
        self.prompt = 'RBCP>'

    def help_load(self):
        print 'load program from file, USAGE: load [file name]'
    def do_load(self, filename):
        for line in fileinput.input(filename):
            line = line.strip()
            self.onecmd(line)

    def do_wrb(self, args):
        global rbcpId
        cmdLine = get_addr_data(args)
        if len(cmdLine) != 2:
            print 'USAGE: wrb [address] [data]'
            return
        else:
            rbcpAddr = hex_int(cmdLine[0])
            rbcpWd = hex_int(cmdLine[1])

            rbcpLen = 1
            rbcpPack = '>B'
            print '- Write %s to 0x%x (ID=%d)' % (rbcpWd,rbcpAddr,rbcpId)
            data = struct.pack(rbcpPack, rbcpWd)
            sitcprbcp.write_registers(ipAddress, rbcpAddr, rbcpLen, data, rbcpId)
            rbcpId += 1
            rbcpId = rbcpId & 255
    def help_wrb(self):
        print 'Write a byte(8bit), USAGE: wrb [address] [data]'

    def do_wrs(self, args):
        global rbcpId
        cmdLine = get_addr_data(args)
        if len(cmdLine) != 2:
            print 'USAGE: wrs [address] [data]'
            return
        else:
            rbcpAddr = hex_int(cmdLine[0])
            rbcpWd = hex_int(cmdLine[1])
            rbcpLen = 2
            rbcpPack = '>H'
            print '- Write %s to 0x%x (ID=%d)' % (rbcpWd,rbcpAddr,rbcpId)
            data = struct.pack(rbcpPack, rbcpWd)
            sitcprbcp.write_registers(ipAddress, rbcpAddr, rbcpLen, data, rbcpId)
            rbcpId += 1
            rbcpId = rbcpId & 255
    def help_wrs(self):
        print 'Write a short(16bit), USAGE: wrs [address] [data]'

    def do_wrl(self, args):
        global rbcpId
        cmdLine = get_addr_data(args)
        if len(cmdLine) != 2:
            print 'USAGE: wrl [address] [data]'
            return
        else:
            rbcpAddr = hex_int(cmdLine[0])
            rbcpWd = hex_int(cmdLine[1])
            rbcpLen = 4
            rbcpPack = '>I'
            print '- Write %s to 0x%x (ID=%d)' % (rbcpWd,rbcpAddr,rbcpId)
            data = struct.pack(rbcpPack, rbcpWd)
            sitcprbcp.write_registers(ipAddress, rbcpAddr, rbcpLen, data, rbcpId)
            rbcpId += 1
            rbcpId = rbcpId & 255
    def help_wrl(self):
        print 'Write a long(32bit), USAGE: wrl [address] [data]'

    def do_rd(self, args):
        global rbcpId
        cmdLine = get_addr_data(args)
        if len(cmdLine) == 2:
            rbcpAddr = hex_int(cmdLine[0])
            rbcpLen = hex_int(cmdLine[1])
            data = sitcprbcp.read_registers(ipAddress, rbcpAddr,rbcpLen, rbcpId)
            adrOfst = rbcpAddr % 16
            adrStart = rbcpAddr - adrOfst
            idx=0
            if len(data) == hex_int(cmdLine[1]):
                print 'RBCP>  Succeed in reading. (ID=%d)' % rbcpId
            else:
                print 'RBCP>  Timeout occuered! (ID=%d)' % rbcpId
            for c in data:
                if idx % 16 ==0:
                    if idx==0:
                        print 'RBCP>  %08x:' % adrStart,
                        for s in range(adrOfst):
                            print '  ',
                        idx = adrOfst
                    else:
                        print 
                        print 'RBCP>  %08x:' % (adrStart + idx),
                elif idx % 8 == 0:
                    print '-',
                print '%02x' % (ord(c)),
                idx +=1
            print
            rbcpId += 1
            rbcpId = rbcpId & 255
        else:
            print 'USAGE: rd [address] [length(byte)]'
    def help_rd(self):
        print 'Read block, USAGE: rd [address] [length(byte)]'

    def do_quit(self, args):
        print '\nBye'
        exit()
    def help_quit(self):
        print 'Quit this program, USAGE: quit'

    def do_cd(self, dirName):
        os.chdir(dirName)
        print ' Current: %s' % os.getcwd()
    def help_cd(self):
        print 'Change working directory, USAGE: cd [directory name]'

    def do_pwd(self, args):
        print ' ' + os.getcwd()
    def help_pwd(self):
        print 'Show the current directory, USAGE: pwd'

    def do_ls(self, dirName):
        if dirName == '':
            dirName = '.'
        cmdLine = os.listdir(dirName)
        for line in cmdLine:
            print line
    def help_ls(self):
        print 'Show the file list in the current directory, USAGE: ls'

def get_addr_data(cmdLine):
    cmdLine = cmdLine.strip()
    tempList = cmdLine.split(' ')
    cmdLine = []
    for s in tempList:
        if s != '':
            cmdLine.append(s)
    return cmdLine

def hex_int(str):
    str = str.strip()
    str = str.lower()

    if str[:2] == '0x' or str[0] == 'x':
        return int(str,16)
    else:
        return int(str)

if __name__ == "__main__":
    argv = sys.argv  # コマンドライン引数を格納したリストの取得
    argc = len(argv) # 引数の個数

    if argc == 1:
        ipAddress = DEFAULT_IP
        rbcpPort = DEFAULT_PORT
    elif argc == 3:
        ipAddress = argv[1]
        rbcpPort = argv[2]
    else:
        print 'Usage: %s <IP address> <Port#>\n\n' % argv[0]
        quit()

    p = rbcpCmd()
    p.cmdloop()

