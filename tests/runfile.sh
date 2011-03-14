#!/bin/bash

TEST_FILE=$1
FUNC_NAME=$2
R=R $R_ARGS

if [ -n "$FUNC_NAME" ]; then
echo "running single function ${FUNC_NAME} in ${TEST_FILE}"
$R --slave <<-EOF
	library('RUnit')
	res <- runTestFile('unit/${TEST_FILE}',
		rngKind='default', rngNormalKind='default',
                testFuncRegexp='${FUNC_NAME}',
                useOwnErrorHandler=TRUE)

        printTextProtocol(res, showDetails=TRUE, traceBackCutOff=0)


EOF

else
echo "running all functions in ${TEST_FILE}"
$R --slave <<-EOF
	library('RUnit')
	res <- runTestFile('unit/${TEST_FILE}',
		rngKind='default', rngNormalKind='default')
	printTextProtocol(res, showDetails=TRUE)
EOF

fi
