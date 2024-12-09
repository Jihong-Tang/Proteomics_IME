#!/bin/bash
TIMESTAMP=$(date)
sed -i "s|TIME_UPDATE|Last updated: $TIMESTAMP|g" README.md