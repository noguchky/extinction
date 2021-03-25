#!/usr/bin/env python

"""
                   ATRY  A7
   ______     ______     ______     ______
  |......|   |..4321|   |......|   |..CCCC|
  |......|   |..8765|   |......|   |..DDDD|
=============================================
                 C   ... SCL
                 D   ... SDI
                 0-7 ...  CS

         enable            disable
CSn1 ... 8b'1000_1000 <--> 8b'1000_1001
CSn2 ... 8b'1001_0000 <--> 8b'1001_0001
CSn3 ... 8b'1001_1000 <--> 8b'1001_1001
 :
CSn8 ... 8b'1100_0000 <--> 8b'1100_0001

enable  = (channel + 8) * 8
disable = (channel + 8) * 8 + 1
"""

import math
import sys
import struct
import socket
import sitcprbcp

def setVth(channel, voltage):
    ip_address = '192.168.10.16'
    addresses  = [ 0x6A, 0x6B, 0x69, 0x68, 0x68, 0x69 ]
    values     = [ 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 ]
    length     = 1
    rbcpid     = 100

    enable  = (channel + 8) * 8
    disable = enable + 1

    values[2] = enable
    values[3] = (int(math.floor(voltage) / 0x100) & 0x3F) + 0x10
    values[4] = (int(math.floor(voltage) % 0x100) & 0xFC)
    values[5] = disable

    try:
        for i in range(len(addresses)):
            print('wrb ' + format(addresses[i], '#04x') + ' ' + format(values[i], '#04x'))
            address = addresses[i]
            data = struct.pack('>B', values[i])
            sitcprbcp.write_registers(ip_address, address, length, data, rbcpid, verify = 0)
    except socket.error as e:
        sys.exit(e)
    except Exception as e:
        sys.exit(e)

if __name__ == '__main__':
    args = sys.argv
    if len(args) < 3:
        sys.exit('usage: python setVth.py [channel] [voltage in milli Volt]')

    try:
        channel = int  (args[1])
        voltage = float(args[2])
    except:
        sys.exit('error: invalid argument, not a number')

    if channel < 1 or channel > 8:
        sys.exit('error: invalid channel (must be between 1-8)')

    if voltage < 0:
        print('warning: voltage must be > 0 mV, voltage set to 0 mV')
        voltage = 0

    if voltage > 3000:
        print('warning: voltage must be > 3000 mV, voltage set to 3000 mV')
        voltage = 3000

    setVth(channel, voltage)
