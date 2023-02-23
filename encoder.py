import unittest

class Encoder():
    @classmethod
    def encode(klass, sequence):
        encoding = None
        if isinstance(sequence, str):
            encoding = [klass.encode_base(base) for base in sequence]
        else:
            raise ValueError(f'Unhandled base: {base}')
        return encoding


class IntegerEncoder(Encoder):
    @staticmethod
    def encode_base(base):
        encoding = None
        if isinstance(base, str):
            if base.upper() == 'A':
                encoding = 1
            elif base.upper() == 'C':
                encoding = 2
            elif base.upper() == 'G':
                encoding = 3
            elif base.upper() == 'T':
                encoding = 4
            else:
                raise ValueError(f'Unhandled base: {base}')
        else:
            raise ValueError(f'Unhandled base: {base}')
        return encoding


class BinaryEncoder(Encoder):
    @staticmethod
    def encode_base(base):
        encoding = None
        if isinstance(base, str):
            if base.upper() == 'A':
                encoding = [0, 0, 0, 1]
            elif base.upper() == 'C':
                encoding = [0, 0, 1, 0]
            elif base.upper() == 'G':
                encoding = [0, 1, 0, 0]
            elif base.upper() == 'T':
                encoding = [1, 0, 0, 0]
            else:
                raise ValueError(f'Unhandled base: {base}')
        else:
            raise ValueError(f'Unhandled base: {base}')
        return encoding

class TestBinaryEncoder(unittest.TestCase):
    def setUp(self):
        self.encoder = BinaryEncoder
        self.base_encoder = self.encoder.encode_base
        self.encoded_a = [0, 0, 0, 1]
        self.encoded_c = [0, 0, 1, 0]
        self.encoded_g = [0, 1, 0, 0]
        self.encoded_t = [1, 0, 0, 0]

    def test_encodeA(self):
        a = 'A'
        self.assertEqual(self.base_encoder(a), self.encoded_a)

    def test_encodeC(self):
        c = 'C'
        self.assertEqual(self.base_encoder(c), self.encoded_c)

    def test_encodeG(self):
        g = 'G'
        self.assertEqual(self.base_encoder(g), self.encoded_g)

    def test_encodeT(self):
        t = 'T'
        self.assertEqual(self.base_encoder(t), self.encoded_t)

    def test_encodeRaisesValueErrorOnUnhandledBase(self):
        with self.assertRaises(ValueError):
            self.base_encoder(None)

    def test_encodeSequence(self):
        acgt = 'ACGT'
        encoded_acgt = [self.encoded_a, self.encoded_c, self.encoded_g, self.encoded_t]
        self.assertEqual(self.encoder.encode(acgt), encoded_acgt)


if __name__ == '__main__':
    unittest.main()
