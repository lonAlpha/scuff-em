#!/usr/bin/env python
#-*- coding: utf-8 -*-

#This script is for gmsh

import socket
import struct
import math
import os
import time
import random

BUFSIZE=1024
HOST='127.0.0.1'
PORT=4000
ADDR=(HOST, PORT)
testSocket = False
logPoint = True

while (True): 
    client=socket.socket(socket.AF_INET, socket.SOCK_STREAM) 
    client.connect(ADDR)

    if testSocket: 
        x=(random.random(), random.random(), random.random())
        client.send(struct.pack("ddd", x[0], x[1], x[2]))
        f = client.recv(BUFSIZE)
        print 'send: ', x[0], x[1], x[2]
        print 'recv: ', struct.unpack('d', f)[0]
#        time.sleep(1)

    else: 
        x = struct.unpack("ddd", os.read(0,24))
        if math.isnan(x[0]):
            break
        client.send(struct.pack("ddd", x[0], x[1], x[2]))
        f = client.recv(BUFSIZE)
        value = struct.unpack('d', f)[0]
        os.write(1,struct.pack("d",value))

        if logPoint: 
            log = open('bgm.log', 'a')
            log.write('Point: '+str(x[0])+ ', '+str(x[1])+', '+str(x[2])+'\n')
            log.write('MeshSize: '+str(value)+'\n')
            log.close()
