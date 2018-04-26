import pytest


def pytest_addoption(parser):
    parser.addoption('--pkg-molecule', action='store_true', default=False, help='test molecule package')
