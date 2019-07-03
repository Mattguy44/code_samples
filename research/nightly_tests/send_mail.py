#!/usr/bin/python3.4

import smtplib  
import sys
import socket

fromaddr = 'moris_Tests@' + socket.gethostname()
toaddrs  = ['matthew.j.ryan@colorado.edu',
            'kurt.maute@colorado.edu',
            'keenan.doble@colorado.edu',
            'mathias.schmidt@colorado.edu',
            'lise.noel@colorado.edu'
            ]
# Removed emails: christian.messe@colorado.edu
server = smtplib.SMTP('localhost')  

for ii in range(0, len(toaddrs)):
  thefrom = "From:" + fromaddr
  theto = "To: " + toaddrs[ii]
  subject = "Subject: " + sys.argv[1]
  body = sys.argv[2]

  msg = "\r\n".join([thefrom,theto,subject,"",body])
  
  print("sending the following message:\n")
  print(msg)
  
  server.sendmail(fromaddr, toaddrs[ii], msg)  
  
server.quit()  
