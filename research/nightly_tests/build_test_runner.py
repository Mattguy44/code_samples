#!/usr/bin/python3.4

# subprocess allows running bash scripts
import subprocess
# os allows python to use environment variables
import os
# re allows for use of regular expressions
import re
# time allows for timekeeping
import time
# Timer executes a command after a period of time on a separate thread
from threading import Timer

# Critical Variables
nightly_test_dir = os.environ['HOME'] + "/Documents/NightlyTests"
package = "moris"
git_repo_url = "ssh://titan/home/git/codes/moris"

# Main Script -------------------------------------------------------------

# Clone the repository
subprocess.call([nightly_test_dir + "/clone_and_compile.sh", nightly_test_dir, package, git_repo_url])

# List build configurations
# MORIS_HAVE_PARALLEL=OFF is not tested because some packages are not yet written with this option in mind. (Oct. 29, 2018)
test_options = ["",
                #"-DMORIS_USE_ARMA=OFF -DMORIS_USE_EIGEN=ON",
                "-DMORIS_HAVE_DEBUG=ON"
                #"-DMORIS_USE_ARMA=OFF -DMORIS_USE_EIGEN=ON -DMORIS_HAVE_DEBUG=ON"
                #"-DMORIS_HAVE_PARALLEL=OFF",
                #"-DMORIS_HAVE_PARALLEL=OFF -DMORIS_HAVE_DEBUG=ON",
                #"-DMORIS_HAVE_PARALLEL=OFF -DMORIS_HAVE_DEBUG=ON -DMORIS_USE_ARMA=OFF -DMORIS_USE_EIGEN=ON"
                ]

test_packages = ["",
                 "-DBUILD_ALG=ON",
                 "-DBUILD_ASR=ON",
                 "-DBUILD_CHR=ON",
                 "-DBUILD_COM=ON",
                 "-DBUILD_CON=ON",
                 "-DBUILD_DLA=ON",
                 "-DBUILD_EXC=ON",
                 "-DBUILD_FEM=ON",
                 "-DBUILD_GEN=ON",
                 "-DBUILD_HMR=ON",
                 "-DBUILD_INT=ON",
                 "-DBUILD_IOS=ON",
                 "-DBUILD_LINALG=ON",
                 "-DBUILD_MAP=ON",
                 "-DBUILD_MDL=ON",
                 #"-DBUILD_MOD=ON",
                 "-DBUILD_MSI=ON",
                 "-DBUILD_MTK=ON",
                 "-DBUILD_NLA=ON",
                 "-DBUILD_OPT=ON",
                 "-DBUILD_SDF=ON",
                 #"-DBUILD_STK=ON",
                 "-DBUILD_TIN=ON",
                 "-DBUILD_TOL=ON",
                 "-DBUILD_XTK=ON",
                 "-DBUILD_ALL=ON"
                 ]

# Run through build tests -------------------------------------------------

# Declare flags and variables across all tests
s_failed_build = 0
s_failed_tests = 0
summary_body = ""

#0#
# Loop through moris packages
for packs in test_packages:
    #1#
    # Reset flags
    failed_build = 0
    failed_tests = 0
    
    # Get package name
    subpackage = re.match("-DBUILD_(.*)=ON", packs)
    if subpackage:
        subpack = subpackage.group(1)
    else:
        subpack = packs
    
    # Set beginnings of error messages
    err_subject = package + " failed to build target " + subpack
    err_summary = "Error summary:"
    err_report = "\n\nError messages:"
    test_err_subject = "Tests under the target " + subpack + " failed."
    test_err_summary = "Test error summary:\n========================================\n\n"
    test_err_report = "\n\nTest error messages:\n========================================\n\n"
    
    # Edge case
    if subpack == "":
        err_subject = package + " failed to build default targets"
    
    # Loop through build options
    for opts in test_options:
        #2#
        # Build the package
        script_output = subprocess.check_output([nightly_test_dir + "/build_package.sh", nightly_test_dir, package, opts, packs], universal_newlines=True)
        
        # Separate CMake output from custom messages
        cmake_status, return_message = script_output.split("@@~~@@")
        print("\n-----\n" + return_message + "\n-----\n")
        
        # Add build status to summary message
        summary_body = summary_body + return_message
        
        # Check if build succeeded
        if "encountered errors" in return_message:
            #3#
            print("FAILED TO BUILD")
            
            # Set flags
            failed_build = 1
            s_failed_build = 1
            
            # Add information to error messages
            log_dir = nightly_test_dir + "/" + package + "/logs"
            err_summary = err_summary + "\n========================================\n\n" + return_message + "\n\n========================================\n"
            err_file = open(log_dir + "/err_output.txt", "r")
            err_report = err_report + "\n========================================\n" + return_message + "\n-----\n\n" + err_file.read()
        else:
            #3#
            # Run tests if build succeeded
            
            # Set timeout parameters
            #timeout = 900 # seconds
            interval = 0.1 # seconds
            
            # Get test process
            test_output = subprocess.check_output([nightly_test_dir + "/run_test.sh", nightly_test_dir, package], universal_newlines=True, timeout=900)
            
            # Send test to timeout manager
            #test_output = timeout_check(test_process, timeout, interval)
            
            # Check if test completed
            if test_output == -9:
                #4#
                print("TEST TIMEOUT")
                test_output = "Test was killed after " + str(timeout) + " seconds. Error code " + str(test_output) + ".\n"
                
                # Set flags
                failed_tests = 1
                s_failed_tests = 1
                
                # Set error message
                test_err_message = "Test timed out with the following build options:\n" + opts + packs + "\n\nTimeout was: " + str(timeout) + " seconds."
                
                # Note status in summary
                summary_body = summary_body + "\nHowever one or more tests timed out.\n"
            else:
                #4#
                # If test completed running, check for failed tests
                if "FAILED" in test_output:
                    #5#
                    print("TEST FAILURE")
                    
                    # Set flags
                    failed_tests = 1
                    s_failed_tests = 1
                    
                    # Set error message
                    test_err_message = "Tests failed with the following build options:\n" +  opts + packs
                    
                    # Note status in summary
                    summary_body = summary_body + "\nHowever one or more tests failed.\n"
                else:
                    #5#
                    # If no failures were found
                    print("TEST SUCCESS")
                    
                    # Note status in summary
                    summary_body = summary_body + "\nAnd all tests were successful.\n"
            
            #3#
            # Set error messages if one or more tests failed
            if failed_tests:
                #4#
                test_err_summary = test_err_summary + test_err_message + "\n\n----------------------------------------\n\n"
                test_err_report = test_err_report + "\n\n----------------------------------------\n\n" + test_err_message + "\n\n-----\n\n" + str(test_output)
        
        #2#
        # Seperate status of each build
        summary_body = summary_body + "\n-----\n"
    
    #1#
    # Send build error message if the build failed
    if failed_build:
        #2#
        # Body of build error email
        err_body = err_summary + err_report
        
        # err_body cannot exceed the maximum string argument length
        if len(err_body) >= 131070:
            #3#
            print("\n\n!!!!!!!!!!!!!!!\n\nError message too long!\n\n!!!!!!!!!!!!!!!\n\n")
            err_body = err_summary + "\n\n!!!!!!!!!!!!!!!\n\nError message was too long to send. Please run with failed configuration and debug manually."
        
        #2#
        # Send build error email
        subprocess.call([nightly_test_dir + "/send_mail.py", err_subject, err_body])
    
    #1#
    # Send test error message if any tests failed
    if failed_tests:
        #2#
        # Body of test error email
        test_err_body = test_err_summary + test_err_report
        
        # test_err_body cannot exceed the maximum string argument length
        if len(test_err_body) >= 131070:
            #3#
            print("\n\n!!!!!!!!!!!!!!!\n\nError message too long!\n\n!!!!!!!!!!!!!!!\n\n")
            test_err_body = test_err_summary + "\n\n!!!!!!!!!!!!!!!\n\nError message was too long to send. Please run with failed configuration and debug manually."
       
        #2#
        # Send test failure email
        subprocess.call([nightly_test_dir + "/send_mail.py", test_err_subject, test_err_body])
#0#

# Summary message ------------------------------------------------

if s_failed_build:
    # Email subject if any builds failed
    summary_subject = "Summary: " + package + " built successfully for some configurations"
else:
    # Email subject if all builds succeeded
    summary_subject = package + " built successfully for all configurations"

if s_failed_tests:
    # Email subject if any tests failed
    summary_subject = summary_subject + ", but did not pass all tests."
else:
    # Email subject if all tests succeeded
    summary_subject = summary_subject + " and passed all tests."

# Send summary email
subprocess.call([nightly_test_dir + "/send_mail.py", summary_subject, summary_body])

print("+-+-+-+-+-+-+-+-+-+ Done +-+-+-+-+-+-+-+-+-+")
