#! /bin/sh
# check script for Arg_parser - POSIX/GNU command line argument parser.
# Copyright (C) 2011-2015 Antonio Diaz Diaz.
#
# This script is free software: you have unlimited permission
# to copy, distribute and modify it.

LC_ALL=C
export LC_ALL
objdir=`pwd`
testdir=`cd "$1" ; pwd`
PARSER="${objdir}"/arg_parser
CPARSER="${objdir}"/carg_parser
framework_failure() { echo "failure in testing framework" ; exit 1 ; }

if [ ! -f "${PARSER}" ] || [ ! -x "${PARSER}" ] ; then
	echo "${PARSER}: cannot execute"
	exit 1
fi
if [ ! -f "${CPARSER}" ] || [ ! -x "${CPARSER}" ] ; then
	echo "${CPARSER}: cannot execute"
	exit 1
fi

if [ -d tmp ] ; then rm -rf tmp ; fi
mkdir tmp
cd "${objdir}"/tmp

in="${testdir}"/test.txt
fail=0

printf "testing arg_parser-%s..." "$2"

for i in "${PARSER}" "${CPARSER}" ; do
	"$i" -a -b 5 -c -carg -o file --orphan file1 file2 > out || fail=1
	cmp ${in} out || fail=1
	printf .

	"$i" -ab5 -c -carg -ofile file1 --orphan file2 > out || fail=1
	cmp ${in} out || fail=1
	printf .

	"$i" --append file1 --block 5 --casual file2 --casual=arg -o file --orphan > out || fail=1
	cmp ${in} out || fail=1
	printf .

	"$i" --append --block=5 --casual= file1 --casual=arg -ofile --orphan file2 > out || fail=1
	cmp ${in} out || fail=1
	printf .

	"$i" -aARG 2> /dev/null
	if [ $? = 1 ] ; then printf . ; else printf - ; fail=1 ; fi

	"$i" -b 2> /dev/null
	if [ $? = 1 ] ; then printf . ; else printf - ; fail=1 ; fi

	"$i" --block 2> /dev/null
	if [ $? = 1 ] ; then printf . ; else printf - ; fail=1 ; fi

	"$i" --block= 2> /dev/null
	if [ $? = 1 ] ; then printf . ; else printf - ; fail=1 ; fi

	"$i" -x 2> /dev/null
	if [ $? = 1 ] ; then printf . ; else printf - ; fail=1 ; fi

	"$i" -u 2> /dev/null
	if [ $? = 3 ] ; then printf . ; else printf - ; fail=1 ; fi
done

echo
if [ ${fail} = 0 ] ; then
	echo "tests completed successfully."
	cd "${objdir}" && rm -r tmp
else
	echo "tests failed."
fi
exit ${fail}
