#!/usr/bin/env python

import h5py
import numpy as np


def load_field(filepath, real, imag=None):
    with h5py.File(filepath, mode='r') as f:
        if f[real].ndim == 2:
            return StaticFieldData(filepath, real, imag)
        elif f[real].ndim == 3:
            return DynamicFieldData(filepath, real, imag)
        else:
            raise ValueError(f'Invalid data "{real}" in {filepath}')


class FieldData:

    def __init__(self, filepath, real, imag=None):
        self.filepath = filepath
        self.real = real
        self.imag = imag

    def open(self):
        self.f = h5py.File(self.filepath, mode='r')

    def close(self):
        del self.f

    def check_data(self):
        if self.real not in self.f:
            raise KeyError(self.real)
        if self.imag is not None:
            if self.imag not in self.f:
                raise KeyError(self.imag)
            rshape = self.f[self.real].shape
            ishape = self.f[self.imag].shape
            if rshape != ishape:
                raise ValueError('Different shape for real and imag')


class StaticFieldData(FieldData):

    def __init__(self, filepath, real, imag=None):
        FieldData.__init__(self, filepath, real, imag)

    def check_data(self):
        ndim = self.f[self.real].ndim
        if ndim != 2:
            raise ValueError(f'Wrong ndim {ndim} (expected 2)')

    @property
    def data(self):
        real = np.asarray(self.f[self.real])
        if self.imag is None:
            return real
        else:
            imag = np.asarray(self.f[self.imag])
            return real + np.complex('j') * imag

    @property
    def shape(self):
        return self.f[self.real].shape


class DynamicFieldData(FieldData):

    def __init__(self, filepath, real, imag=None):
        FieldData.__init__(self, filepath, real, imag)
        self.i = 0

    def check_data(self):
        ndim = self.f[self.real].ndim
        if ndim != 3:
            raise ValueError(f'Wrong ndim {ndim} (expected 3)')

        ndim = self.f[self.real].ndim
        if ndim != 3:
            raise ValueError(f'Wrong ndim {ndim} (expected 3)')

    @property
    def shape(self):
        return self.f[self.real].shape[:2]

    def __len__(self):
        return self.f[self.real].shape[2]

    def __getitem__(self, n):
        real = np.asarray(self.f[self.real][:, :, n])
        if self.imag is None:
            return real
        else:
            imag = np.asarray(self.f[self.imag][:, :, n])
            return real + np.complex('j') * imag

    def get_min(self, value='real'):
        if value == 'real':
            return np.min([np.min(v.real) for v in self])
        elif value == 'imag':
            return np.min([np.min(v.imag) for v in self])
        elif value == 'angle':
            return np.min([np.min(np.angle(v)) for v in self])
        elif value == 'abs':
            return np.min([np.min(np.abs(v)) for v in self])
        else:
            raise KeyError(value)

    def get_max(self, value='real'):
        if value == 'real':
            return np.max([np.max(v.real) for v in self])
        elif value == 'imag':
            return np.max([np.max(v.imag) for v in self])
        elif value == 'angle':
            return np.max([np.max(np.angle(v)) for v in self])
        elif value == 'abs':
            return np.max([np.max(np.abs(v)) for v in self])
        else:
            raise KeyError(value)

    def __iter__(self):
        self.i = 0
        return self

    def __next__(self):
        if self.i < len(self):
            v = self[self.i]
            self.i += 1
            return v
        else:
            raise StopIteration
