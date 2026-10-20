"""Dump the chunk directory of a SeqBuilder .sbd file (research aid)."""
import struct, sys

def chunks(d):
    objtab = struct.unpack('>I', d[0x2c:0x30])[0]
    out = []
    o = objtab
    while o + 0x40 <= len(d):
        idx, _, typ, sub, _, grp = struct.unpack('>II4s4sII', d[o:o + 24])
        name = d[o + 25:o + 25 + d[o + 24]].decode('latin1')
        cc, n, _, _, first = struct.unpack('>5I', d[grp:grp + 20])
        for k in range(n):
            e = first + 32 * k
            t, cid, fl, ln, _, _, _, off = struct.unpack('>4s7I', d[e:e + 32])
            out.append((typ.decode() + '/' + name, t.decode('latin1'), cid, off, ln))
        o += 0x40
    return out

if __name__ == '__main__':
    d = open(sys.argv[1], 'rb').read()
    for c in chunks(d):
        print(*c)

def props(b, depth=0, out=None, path=''):
    """CommonBlockFmt property tree: 8-byte key, u32 LE length, value (nested when it parses as props)."""
    out = [] if out is None else out
    i = 0
    if b.startswith(b'CommonBlockFmt\0'):
        i = 17
    while i + 12 <= len(b):
        key = b[i:i + 8].rstrip(b'\0').decode('latin1')
        ln = struct.unpack('<I', b[i + 8:i + 12])[0]
        val = b[i + 12:i + 12 + ln]
        if ln > len(b) - i - 12 or not key.isprintable():
            out.append((path + '?', b[i:].hex()[:80])); break
        out.append((path + key, ln, val[:24].hex()))
        if ln >= 12 and val[:1].isalpha() and val[8:12] and struct.unpack('<I', val[8:12])[0] <= ln - 12:
            props(val, depth + 1, out, path + key + '.')
        i += 12 + ln
    return out
