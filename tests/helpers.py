#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy


def array_1d_test(function, input, reference, rel_error, add_msg=""):
    """helper function to test a function with one input parameter"""
    res = numpy.zeros(len(reference))

    for i, in1 in enumerate(input):
        res[i] = function(in1)

    error = numpy.sum(numpy.absolute((res - reference) / reference))

    msg = f"Array test of {function.__name__} failed with error {error}, allowed is {rel_error}."
    if add_msg:
        msg += f" {add_msg}"

    assert error < rel_error, msg


def array_2d_test(function, inputs, reference, rel_error, add_msg=""):
    """helper function to test a function with two input parameter"""
    res = numpy.zeros(len(reference))

    for i, (in1, in2) in enumerate(zip(inputs[0], inputs[1])):
        res[i] = function(in1, in2)

    error = numpy.sum(numpy.absolute((res - reference) / reference))

    msg = f"Array test of {function.__name__} failed with error {error}, allowed is {rel_error}."
    if add_msg:
        msg += f" {add_msg}"

    assert error < rel_error, msg
