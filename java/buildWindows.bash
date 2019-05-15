#!/bin/bash

pushd $(dirname "$0") > /dev/null

export CLASSPATH=`find \lib -name "*.jar" | tr "\n" ";"`

find . -name "*.java" | xargs javac -Xlint:unchecked

CLASS=`find . -name "*.class" | tr "\n" " "`

jar cf astria.jar $CLASS

popd > /dev/null

mv astria.jar ../lib
cd org/astria
find . -name "*.class" | xargs rm