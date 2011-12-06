#!/usr/bin/env python
# 
# NAME
#
#	dgmsocket class
#
# DESCRIPTION
#
#	The 'dgmsocket' class provides a very simple wrapper the standard
#	python socket API.
#
#	More specifically, this class provides datagram socket services.
#
# HISTORY
#
# 25 March 2006
# o Initial development implementation
#
# 06 December 2011
# o Clean-up the socket communication
#
import socket

class C_dgmsocket :
	# 
	# Member variables
	#
	# 	- Core variables
	mstr_obj	= 'C_dgmsocket';	# name of object class
        mstr_name	= 'void';		# name of object variable
        m_id		= -1; 			# id of agent
        m_iter		= 0;			# current iteration in an
                                		# 	arbitrary processing 
						#	scheme
        m_verbosity	= 0;			# debug related value for 
						#	object
        m_warnings	= 0;              	# show warnings 
						#	(and warnings level)
	
	#
	#	- Class variables
	m_dgmsocket	= None;
	mstr_remoteHost	= 'localhost';
	m_port		= 1701;
	
	#
	# Methods
	#
	# Core methods - construct, initialise, id
	def core_construct(	self,
				astr_obj	= 'C_dgmsocket',
				astr_name	= 'void',
				a_id		= -1,
				a_iter		= 0,
				a_verbosity	= 0,
				a_warnings	= 0) :
		self.mstr_obj		= astr_obj
		self.mstr_name		= astr_name
		self.m_id		= a_id
		self.m_iter		= a_iter
		self.m_verbosity	= a_verbosity
		self.m_warnings		= a_warnings
	def __str__(self):
		print 'mstr_obj\t\t= %s' 	% self.mstr_obj
		print 'mstr_name\t\t= %s' 	% self.mstr_name
		print 'm_id\t\t\t= %d' 		% self.m_id
		print 'm_iter\t\t\t= %d'	% self.m_iter
		print 'm_verbosity\t\t= %d'	% self.m_verbosity
		print 'm_warnings\t\t= %d'	% self.m_warnings
		return 'This class provides a *very* simple wrapper framework about datagram sockets.'
	def __init__(self, astr_hostname = 'localhost', a_port = 1701):
		self.core_construct()
		self.mstr_remoteHost	= astr_hostname
		self.m_port 		= a_port
		self.m_dgmsocket = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
	def tx(self, str_payload):
		self.m_dgmsocket.sendto(str_payload, (self.mstr_remoteHost, self.m_port))

