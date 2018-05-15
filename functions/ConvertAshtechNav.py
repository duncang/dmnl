# (c) Damien Dusha 2006
# Queensland University of Technology
# All Rights Reserved
#
# CVS Information:
# $Id: ConvertAshtechNav.py 807 2007-04-12 23:03:11Z greerd $
#
# $Log$
#

# Import the oeprating system libraries
import os
import shutil
import sys
import getopt
import time
import calendar
import locale

#----------------------------------------------------------------------------------------
# Global Variables (default for many things)
InputRoot = "./"
OutputRoot = "./Images/"
OutputFileName = "OutputFile.txt"
UseOutputFile = 0
Extension = ".raw"



#----------------------------------------------------------------------------------------
def RecursiveFullPathList(CurrentDirectory, FileList=[]):

	# Strip whitespace off the argument
	CurrentDirectory = CurrentDirectory.strip()
	
	# Strip off any trailing slash
	CurrentDirectory = CurrentDirectory.rstrip("/")
	CurrentDirectory = CurrentDirectory.rstrip("\\")
	
	# Add on a trailing slash
	CurrentDirectory = CurrentDirectory + "/"
	
	# List the current directory
	Files = os.listdir(CurrentDirectory)
	
	# Sort the files
	Files.sort()
	
	# Go through each file in the list....
	for Filename in Files:
		# Create the full path for this file
		FullPath = CurrentDirectory + Filename
		
		# Check if the file is a directory
		if os.path.isdir(FullPath):
			
			#Recursively list all the files.....
			RecursiveFullPathList(FullPath, FileList)
			
		else:
		
			# Write the full path to the stdout
			FileList.append(FullPath)
			
	return FileList


#----------------------------------------------------------------------------------------
def SplitImageFileName(FileName):

	# Split the directory and File Name Parts
	BreakupFile = FileName.rsplit("/", 1)
	
	# Assign each of the parts
	File = BreakupFile[1]
	Directory = BreakupFile[0] + "/"
	
	# Split the file into its body and file Extension
	BreakupFile = File.rsplit(".",1)
	Body = BreakupFile[0]
	Extension = "." + BreakupFile[1]
	
	# Split the body into its date and time components
	# Two possible forms of the filename
	# Prefix-yyyymmdd-hhmmss-mmm.xyz
	# Prefix-nnnnnnnnnn.raw
	BreakupFile = Body.split("-")
	
	if (len(BreakupFile) < 3):
	
		# It's the raw sequence number
		Prefix = BreakupFile[0]
		Sequence = int(BreakupFile[1])
		
		# Set to an arbitary "DateTime.MinValue" from .NET
		Year = 0001
		Month = 1
		Day = 1
		Hour = 0
		Minute = 0
		Second = 0
		Millisecond = 0
		
	else:
	
		# Full date-time
		Prefix = BreakupFile[0]
		Year = int(BreakupFile[1]) / 10000
		Month = (int(BreakupFile[1]) % 10000) / 100
		Day = (int(BreakupFile[1]) % 100)
		Hour = int(BreakupFile[2]) / 10000
		Minute = (int(BreakupFile[2]) % 10000) / 100
		Second = (int(BreakupFile[2]) % 100)
		Millisecond = int(BreakupFile[3])
		
	# return the unix time
	return [Directory, Body, Extension, Prefix, Year, Month, Day, Hour, \
			Minute, Second, Millisecond, DateTimeToUnixTme(Year, Month, Day, \
			Hour, Minute, Second)]

#----------------------------------------------------------------------------------------		
def DateTimeToUnixTme(Year, Month, Day, Hour, Minute, Second, TimeZone=0):

	# Get the weekday (0 is monday)
	WeekDay = calendar.weekday(Year, Month, Day)
	
	# Get the day of the year as a decimal number
	
	
	# Create the 9-tuple for getting unix time (-1 for system chooses daylight savings
	# Note that the time must be LOCAL TIME
	TimeTuple = (Year, Month, Day, Hour, Minute, Second, WeekDay, 0 , -1)
	
	# Get the unix time out of the tuple
	return int(time.mktime(TimeTuple))

#----------------------------------------------------------------------------------------	
def usage():

	# Print out the usage
	sys.stdout.write("Usage Instructions:\n")
	sys.stdout.write("-------------------\n\n")
	sys.stdout.write(sys.argv[0] + " uses the current directory, and outputs to stdout, searches for \".raw files\"\n")
	sys.stdout.write("-f --filename,   Output the contents to a file\n")
	sys.stdout.write("-d --directory,  Chooses the root directory to search for files\n")
	sys.stdout.write("-e --extension,  Extension to search for\n")
	

#----------------------------------------------------------------------------------------
# Main Script

# Get the arguments to the script
try:                                
    opts, args = getopt.getopt(sys.argv[1:], "f:d:e:", ["filename=", "directory=", "extension="]) 
except getopt.GetoptError:           
	usage()                          
	sys.exit(2)   
	
for opt, arg in opts:

	if opt in ("-f", "--filename"):
		OutputFileName = arg
		UseOutputFile = 1
		
	elif opt in ("-d", "--directory"):
		InputRoot = arg
		
	elif opt in ("-e", "--extension"):
		Extension = arg
		
	else:
		sys.stdout.write("Unrecognised option\n")
		usage()
		system.exit(2)
		


# Get the list of files
FileList = RecursiveFullPathList(InputRoot)

# Sort the list of files
FileList.sort()

# Remove all that do not match the file type
for FileName in FileList[:]:	# make a copy of the list
	
	# Check if it is a "raw" file
	if (FileName.find(Extension) == -1):
		FileList.remove(FileName)
		
# Open the output file, if nessessary
if (UseOutputFile == 1):
	try:
		OutputFile = open(OutputFileName, "w")
	except:
		sys.stdout.write("Error Opening File\n")

	

# Send the list of files out to stdio
for FileName in FileList:
	OutFileName = FileName + ".e"
	os.system("python d_replace.py " + FileName + " " + OutFileName)

