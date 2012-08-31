#include "translate_cigar.h"
#include "bwtaln.h"

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;

namespace {
    const char* OPS = "MIDSN";
    int cigar_len(bwa_cigar_t* cigar, int n) {
        int len = 0;
        for (int i = 0; i < n; ++i) {
            unsigned short op = __cigar_op(cigar[i]);
            if (op <= 1)
                len += __cigar_len(cigar[i]);
        }
        return len;
    }
}

struct CigarBuilder {
    CigarBuilder()
        : cigar(0)
        , len(0)
        , cap(0)
    {
    }

    ~CigarBuilder() {
        if (cigar)
            free(cigar);
    }

    bwa_cigar_t* release() {
        bwa_cigar_t* rv = cigar;
        cigar = 0;
        return rv;
    }

    void push(char op, int len) {
        push(__cigar_create(op, len));
    }

    void push(bwa_cigar_t cig) {
        if (len && __cigar_op(cig) == __cigar_op(cigar[len-1])) {
            int newlen = __cigar_len(cig)+__cigar_len(cigar[len-1]);
            cigar[len-1] = __cigar_create(__cigar_op(cig), newlen);
        } else {
            if (cap == len) {
                ++cap;
                cap *= 2;
                cigar = (bwa_cigar_t*)realloc(cigar, cap*sizeof(bwa_cigar_t));
                if (cigar == NULL) {
                    cerr << "CigarBuilder: out of memory.\n";
                    exit(1);
                }
            }
            cigar[len++] = cig;
        }
    }

    bwa_cigar_t* cigar;
    int len;
    int cap;
};


struct CigarTranslator {
    CigarTranslator(
            const char* seq_cigar,
            uint32_t start_pos, // start pos on seq
            bwa_cigar_t* read_cigar,
            int n_cigar,
            int total_read_len
            )
        : seq_cigar(seq_cigar)
        , p_seq_cigar(seq_cigar)
        , start_pos(start_pos)
        , read_cigar(read_cigar)
        , read_cigar_idx(0)
        , n_cigar(n_cigar)
        , total_read_len(total_read_len)
        , cpos(0)
    {
        seq_advance();
        read_advance();
    }

    int tr_seqop(char op) {
        switch (op) {
            case 'M': return 0; break;
            case 'I': return 1; break;
            case 'D': return 2; break;
            case 'S': return 3; break;
            case 'N': return 4; break;
            default:
                throw runtime_error(string("Unknown cigar operation: ") + op);
                break;
        }
    }

    void in_match() {
/*
        if (seq_len >= read_len) {
            cb.push(read_op, read_len);
            seq_len -= read_len;
            read_len = 0;
        } else {
            cb.push(read_op, seq_len);
            read_len -= seq_len;
            seq_len = 0;
        }
*/

        switch (OPS[read_op]) {
            case 'M':
            case 'N':
            case 'D':
                if (seq_len >= read_len) {
                    cb.push(read_op, read_len);
                    seq_len -= read_len;
                    read_len = 0;
                } else {
                    cb.push(read_op, seq_len);
                    read_len -= seq_len;
                    seq_len = 0;
                }
                break;

            case 'I':
                cb.push(read_op, read_len);
                read_len = 0;
                break;

            default:
                throw runtime_error("Unknown cigar op in read");
                break;
        }
    }

    void in_insertion() {
        switch (OPS[read_op]) {
            case 'M':
                if (seq_len < read_len) {
                    cb.push(1, seq_len);
                    read_len -= seq_len;
                    seq_len = 0;
                } else {
                    cb.push(1, read_len);
                    seq_len -= read_len;
                    read_len = 0;
                }
                break;

            case 'I':
                cb.push(read_op, read_len);
                read_len = 0;
                break;

            case 'N':
            case 'D':
                
                if (seq_len > read_len) {
                    seq_len -= read_len;
                    read_len = 0;
                } else {
                    read_len -= seq_len;
                    seq_len = 0;
                }
                break;

            default:
                throw runtime_error("Unknown cigar op in read");
                break;
        }
    }

    void in_deletion() {
        switch (OPS[read_op]) {
            case 'M':
                cb.push(tr_seqop(seq_op), seq_len);
                seq_advance();
                break;

            case 'I':
                cb.push(tr_seqop(seq_op), seq_len);
                seq_advance();

                cb.push(read_op, read_len);
                read_advance();
                break;

            case 'N':
            case 'D':
                cb.push(tr_seqop(seq_op), seq_len);
                seq_len = 0;
                break;
            default:
                throw runtime_error("Unknown cigar op in read");
                break;
        }
    }

    void exec() {
        find_start_pos();
        if (read_cigar == NULL) {
            int len = 0;
            while (len < total_read_len && !eos()) {
                int dist = total_read_len - len;
                if (seq_len < dist) {
                    cb.push(tr_seqop(seq_op), seq_len);
                    len += seq_len;
                    seq_advance();
                } else {
                    cb.push(tr_seqop(seq_op), dist);
                    break;
                }
            }
            return;
        }

        while (!eor() && !eos()) {
            if (seq_len == 0) seq_advance();
            if (read_len == 0) read_advance();
            if (OPS[read_op] == 'S') {
                cb.push(read_op, read_len);
                read_len = 0;
                if (!eor())
                    read_advance();
                continue;
            }

            switch (seq_op) {
                case '=':
                case 'M':
                case 'X':
                    in_match();
                    break;

                case 'I':
                    in_insertion();
                    break;

                case 'N':
                case 'D':
                    in_deletion();
                    break;

                default:
                    throw runtime_error(string("Invalid cigar character: ") + seq_op);
            }
        }

        while (!eor()) {
            if (read_len == 0)
                read_advance();

            if (OPS[read_op] == 'M' || OPS[read_op] == 'I' || OPS[read_op] == 'S')
                cb.push(tr_seqop('S'), read_len);
            read_len = 0;
        }
    }

    void find_start_pos() {
        while (cpos < start_pos && !eos()) {
            if (seq_len == 0) seq_advance();
            int dist = start_pos - cpos;
            switch (seq_op) {
                case '=':
                case 'M':
                case 'X':
                case 'I':
                    if (seq_len > dist) {
                        seq_len -= start_pos - cpos;
                        cpos = start_pos;
                    } else {
                        cpos += seq_len;
                        seq_len = 0;
                    }
                    break;

                case 'N':
                case 'D':
                    seq_len = 0;
                    break;

                default:
                    throw runtime_error(string("Invalid cigar character: ") + seq_op);
                    break;
            }
        }
        if (cpos < start_pos) {
            stringstream ss;
            ss << "Failed to seek to position " << start_pos << " in cigar string '" << seq_cigar << "'";
            throw runtime_error(ss.str());
        }
    }

    bool eos() const {
        return seq_len == 0 && *p_seq_cigar == 0;
    }

    bool eor() const {
        return read_len == 0 && read_cigar_idx >= n_cigar;
    }

    void seq_advance() {
        char* end;
        seq_len = strtoul(p_seq_cigar, &end, 10);
        p_seq_cigar = end;
        seq_op = *p_seq_cigar++;
    }

    void read_advance() {
        if (!read_cigar)
            return;

        read_len = __cigar_len(read_cigar[read_cigar_idx]);
        read_op = __cigar_op(read_cigar[read_cigar_idx++]);
    } 

    
    const char* seq_cigar;
    const char* p_seq_cigar;
    uint32_t start_pos;
    bwa_cigar_t* read_cigar;
    int read_cigar_idx;
    int n_cigar;
    int total_read_len;
    uint32_t cpos;

    char seq_op;
    int seq_len;
    int read_op;
    int read_len;
    CigarBuilder cb;
};

bwa_cigar_t* translate_cigar(
    const char* cigar,
    uint32_t start,
    bwa_cigar_t* read_cigar,
    int n_cigar,
    int read_len,
    int* n_cigar_out)
{
    CigarTranslator ct(cigar, start, read_cigar, n_cigar, read_len);
    try {
        ct.exec();
        *n_cigar_out = ct.cb.len;
        return ct.cb.release();
    } catch (const exception& e) {
        cerr << "Error translating cigar string: " << e.what() << "\n";
        return NULL;
    }
}

#if 0
void print_cigar(bwa_cigar_t* cigar, int len) {
    for (int i = 0; i < len; ++i) {
        cout << __cigar_len(cigar[i]) << OPS[__cigar_op(cigar[i])];
    }
    cout << "\n";
}

int main(int argc, char** argv) {
    if (argc != 4) {
        cerr << "Give 2 cigar strings, pos\n";
        return 1;
    }

    uint32_t pos = atoi(argv[3]);
    const char* seq_cigar = argv[1];
    CigarBuilder cb;
    const char* p = argv[2];
    static char trans[256] = {0};
    trans['M'] = 0;
    trans['I'] = 1;
    trans['D'] = 2;
    trans['S'] = 3;
    while (*p) {
        char* end;
        int len = strtoul(p, &end, 10);
        p = end;
        cb.push(trans[*p++], len);
    }

    print_cigar(cb.cigar, cb.len);

    int retlen;
    int len = cigar_len(cb.cigar, cb.len);
    bwa_cigar_t* ret = translate_cigar(argv[1], pos, cb.cigar, cb.len, len, &retlen);
    if (ret) {
        cout << "RESULT: ";
        print_cigar(ret, retlen);
    } else {
        cerr << "NO RET!\n";
    }

    ret = translate_cigar(argv[1], pos, NULL, 0, 10, &retlen);
    cout << "RESULT2: ";
    print_cigar(ret, retlen);

    return 0;
}
#endif
