import re
from collections import Counter

####################
# Functions


def cs_to_list(cs_string):
    '''
    Convert the cs_strings into a list with values and genomic positions (0-based).

    Parameters:
        cs_string - string, the cs tag string from BAM file mapped with minimap2 with --cs.

    Return: list of separate values in cs string. Item with format: [low, high, cs tag, cs value].
    '''
    # 0 based
    pos = 0
    # length functions dictionary
    cslenfuncs = {
        ':': int,
        '*': lambda x: 1,
        '+': lambda x: 0,
        '-': len,
        '~': lambda x: int(
            re.sub('[a-z]', '', x)
        )
    }
    cs_mark = re.sub(
        '[0-9a-z]', ' ', cs_string
    ).strip().split()
    cs_value = re.sub(
        '[:*\-+~]', ' ', cs_string
    ).strip().split()
    cslist = list()
    for a, b in zip(cs_mark, cs_value):
        low = pos
        pos += cslenfuncs[a](b)
        high = pos
        cslist.append([low, high, a, b])
    return cslist


def cigar_to_list(cigar_string):
    '''
    Convert the cigar_strings into a list with values and genomic positions (0-based).

    Parameters:
        cs_string - string, the cs tag string from BAM file mapped with minimap2 with --cs.

    Return: list of separate values in cs string. Item with format: [low, high, cigar tag, cigar value]
    '''
    # 0 based
    pos = 0
    cigarlenfuncs = {
        'M': int,
        'D': int,
        'N': int,
        'S': lambda x: 0,
        'I': lambda x: 0
    }
    cigar_mark = re.sub(
        '[0-9]+', ' ', cigar_string
    ).strip().split()
    cigar_num = re.sub(
        '[a-zA-Z]', ' ', cigar_string
    ).strip().split()
    cigarlist = list()
    for a, b in zip(cigar_mark, cigar_num):
        low = pos
        pos += cigarlenfuncs[a](b)
        high = pos
        cigarlist.append([low, high, a, b])
    return cigarlist


def md_to_list(md_string):
    digit_char = [
        '0', '1', '2', '3', '4',
        '5', '6', '7', '8', '9'
    ]
    up_char = [
        'A', 'B', 'C', 'D', 'E', 'F', 'G',
        'H', 'I', 'J', 'K', 'L', 'M', 'N',
        'O', 'P', 'Q', 'R', 'S', 'T',
        'U', 'V', 'W', 'X', 'Y', 'Z'
    ]
    md_split = re.findall(
        '([0-9]+|[A-Z]|\^[A-Z]+)', md_string
    )
    mdlist = list()
    for a in md_split:
        if a[0] == '^':
            mdlist.append(
                [a, 'D', len(a) - 1]
            )
        elif a[0] in digit_char:
            mdlist.append(
                [a, 'M', int(a)]
            )
        elif a[0] in up_char:
            mdlist.append(
                [a, 'X', len(a)]
            )
        else:
            pass
    return mdlist


def merge_cigar_md(cigar_string,
                   md_string):
    # cigar ignores mismatch.
    # md ignores insertion and splicing, and the matched region
    # flanking insertion and splicing will merge togather
    # merged list will remove soft clip
    cigar = cigar_to_list(cigar_string)
    # cigar with items: ref_low, ref_high, cigar_mark, cigar_value
    # remove soft clip
    cigar = [a for a in cigar if a[2] != 'S']
    md = md_to_list(md_string)
    # md with items: md_value, md_type, length

    ####################
    # position of read sequence by cigar and md
    # sequence position by cigar includeing M, I, (S), but not D
    # sequence position by md including M, X, but not D, I and S

    # cigar
    # 0 based
    current_cigar_length = 0
    current_md_length = 0
    current_cigar_start = 0
    current_md_start = 0
    # add more columns to cigar:
    #   ref_low, ref_high, cigar_mark, cigar_value,
    #   cigar_low, cigar_high, md_low, md_high
    for i, item in enumerate(cigar):
        low, high, cigar_mark, cigar_value = item
        if cigar_mark == 'M':
            current_cigar_length += int(cigar_value)
            current_md_length += int(cigar_value)
        elif cigar_mark == 'I':
            current_cigar_length += int(cigar_value)
        else:
            pass
        # use 0-base
        cigar[i] = item + [
            current_cigar_start, current_cigar_length,
            current_md_start, current_md_length
        ]
        current_cigar_start = current_cigar_length
        current_md_start = current_md_length
    # md
    # 0 based
    current_md_length = 0
    current_md_start = 0
    # add more columns to md:
    #   md_value, md_type, length,
    #   md_low, md_high
    for i, item in enumerate(md):
        md_value, md_type, md_len = item
        if md_type in ['M', 'X']:
            current_md_length += md_len
        else:
            pass
        # use 0-base
        md[i] = item + [current_md_start, current_md_length]
        current_md_start = current_md_length

    # find mismatches in md
    mismatches = [
        md[i] for i in range(len(md))
        if md[i][1] == 'X'
    ]
    # find deletions in md
    deletions = [
        md[i] for i in range(len(md))
        if md[i][1] == 'D'
    ]
    # add mismatches into CIGAR
    new_cigar = list()
    for i in range(len(cigar)):
        md_start = cigar[i][6]
        md_end = cigar[i][7]
        # get all the mismatches in the region
        inside_mismatches = [
            a for a in mismatches
            if a[3] >= md_start and a[3] < md_end
        ]
        if len(inside_mismatches) == 0:
            if cigar[i][2] == 'D':
                # process deletions
                deletion_item = [
                    a for a in deletions
                    if a[3] == md_start
                ][0]
                new_cigar.append(
                    cigar[i][0:3] +
                    [deletion_item[0].replace('^', '')]
                )
            else:
                new_cigar.append(cigar[i][0:4])
        else:
            ref_start = cigar[i][0]
            ref_end = cigar[i][1]
            # give insite_mismatches more columns:
            # md_value, md_type, length,
            # md_low, md_high, ref_start, ref_end
            for j, item in enumerate(inside_mismatches):
                # mismatch position is 0 based
                # ref position is 0 based
                dist_md = item[3] - md_start
                thepos = dist_md + ref_start
                inside_mismatches[j] = item + [thepos, thepos + 1]
            new_items = list()
            the_start = ref_start
            for item in inside_mismatches:
                item_ref_start = item[5]
                item_ref_end = item[6]
                new_items.append(
                    [item_ref_start, item_ref_end, 'X', item[0]]
                )
                if (item_ref_start - the_start) >= 1:
                    new_items.append(
                        [the_start, item_ref_start,
                         'M', '{}'.format(item_ref_start - the_start)]
                    )
                the_start = item_ref_end
            if (ref_end - the_start) >= 1:
                # the last part
                new_items.append(
                    [the_start, ref_end,
                     'M', '{}'.format(ref_end - the_start)]
                )
            else:
                pass
            new_items = sorted(new_items, key = lambda a: a[0])
            new_cigar.extend(new_items)
    return new_cigar


def cigar_md_to_read_coordinate(cigar_md, reference_start = 0):
    read_coord = list()
    read_pos = 0
    for ref_low, ref_high, cigar_mark, cigar_value in cigar_md:
        ref_len = ref_high - ref_low
        if cigar_mark in ['M', 'X']:
            read_coord.append(
                [read_pos, read_pos + ref_len, cigar_mark, cigar_value]
            )
            read_pos += ref_len
        elif cigar_mark == 'I':
            seq_len = int(cigar_value)
            read_coord.append(
                [read_pos, read_pos + seq_len, cigar_mark, cigar_value]
            )
            read_pos += seq_len
        elif cigar_mark in ['D', 'N']:
            read_coord.append(
                [read_pos, read_pos, cigar_mark, cigar_value]
            )
    return read_coord


def convert_cigar_md_to_cs_list(cigar_string,
                                md_string,
                                read_seq,
                                ref_seq):
    # deal with soft_clip
    cigar = cigar_to_list(cigar_string)
    if cigar[0][2] == 'S':
        # process head soft clip
        soft_len = int(cigar[0][3])
        read_seq = read_seq[soft_len:]
    else:
        pass
    if cigar[-1][2] == 'S':
        # process tail soft clip
        soft_len = int(cigar[-1][3])
        read_seq = read_seq[:-1 * soft_len]
    else:
        pass
    # merge cigar with md
    cm = merge_cigar_md(cigar_string, md_string)
    # convert marks
    mark_cigar_cs_dict = {
        'M': ':', 'X': '*',
        'I': '+', 'D': '-', 'N': '~'
    }
    new_cs = [
        [a[0], a[1],
         mark_cigar_cs_dict[a[2]], '{}'.format(a[3])]
        for a in cm
    ]
    # read_seq position by cigar
    cigar_pos = list()      # sequence position by cigar
    seq_pos = 0
    for low, high, cigar_mark, cigar_value in cm:
        if cigar_mark == 'M':
            seq_pos += int(cigar_value)
        elif cigar_mark == 'X':
            seq_pos += 1
        elif cigar_mark == 'I':
            seq_pos += int(cigar_value)
        elif cigar_mark == 'S':
            seq_pos += int(cigar_value)
        else:
            pass
        # convert to 0-base
        cigar_pos.append(seq_pos - 1)
    ##########
    # deletion to lower
    for i in range(len(new_cs)):
        if new_cs[i][2] == '-':
            old_val = new_cs[i][3]
            new_cs[i][3] = old_val.lower()
        else:
            pass
    ##########
    # add insertion sequence
    for i in range(len(new_cs)):
        if new_cs[i][2] == '+':
            i_end = cigar_pos[i] + 1
            old_val = int(new_cs[i][3])
            i_start = i_end - old_val
            i_seq = read_seq[
                i_start : i_end
            ].lower()
            new_val = i_seq
            new_cs[i][3] = new_val
        else:
            pass
    ##########
    # add mismatch of the read_seq
    for i in range(len(new_cs)):
        if new_cs[i][2] == '*':
            old_val = new_cs[i][3]
            new_val = (
                old_val + read_seq[cigar_pos[i]]
            ).lower()
            new_cs[i][3] = new_val
        else:
            pass
    ##########
    # add splicing pattern
    for i in range(len(new_cs)):
        if new_cs[i][2] == '~':
            intron_start = new_cs[i][0]
            intron_end = new_cs[i][1]
            old_val = new_cs[i][3]
            split1 = ref_seq[
                intron_start : intron_start + 2
            ].lower()
            split2 = ref_seq[
                intron_end - 2 : intron_end
            ].lower()
            new_val = '{}{}{}'.format(
                split1, old_val, split2
            )
            new_cs[i][3] = new_val
        else:
            pass
    return new_cs

####################
# Objects


class CS(object):
    '''
    Store CS tags in the SAM/BAM file generated by minimap2 with --cs parameter.

    Also accept CIGAR, MD, read_seq, ref_seq as input.
    '''
    name = 'CS'

    def __init__(self, cs_tag_list, contig, start_pos, strand):
        super(CS, self).__init__()
        self._cs = cs_tag_list
        self._start_pos = start_pos
        self._contig = contig
        self._strand = strand
        self._read_length = self._get_read_length()
        self._intron_count = self._count_intron()
        self._splice_pair_count = self._count_splice_pair()
        (self._splice_l_count,
         self._splice_r_count) = self._count_splice_site()
        self._mismatch_count = self._count_mismatch()
        self._mismatch_type_count = self._count_mismatch_type()
        self._insertion_count = self._count_insertion()
        self._insertion_length = self._calculate_insertion_length()
        self._deletion_count = self._count_deletion()
        self._deletion_length = self._calculate_deletion_length()

    @classmethod
    def from_cs_tag_list(cls, cs_tag_list,
                         contig, start_pos, strand):
        new_cs_object = CS(
            cs_tag_list,
            contig = contig,
            start_pos = start_pos,
            strand = strand
        )
        return new_cs_object

    @classmethod
    def from_cs_tag_string(cls, cs_tag_string,
                           contig, start_pos, strand):
        cs = cs_to_list(cs_tag_string)
        return cls.from_cs_tag_list(
            cs,
            contig = contig,
            start_pos = start_pos,
            strand = strand
        )

    @classmethod
    def from_cigar_string(cls, cigar_string, md_string,
                          read_seq, ref_seq,
                          contig, start_pos, strand):
        cs = convert_cigar_md_to_cs_list(
            cigar_string, md_string,
            read_seq, ref_seq
        )
        return cls.from_cs_tag_list(
            cs,
            contig = contig,
            start_pos = start_pos,
            strand = strand
        )

    def get_cs_tag_string(self):
        return ''.join(
            ['{}{}'.format(c, d) for a, b, c, d in self._cs]
        )

    def get_relative_position(self):
        return self._cs

    def get_contig_position(self):
        return [
            [a + self._start_pos,
             b + self._start_pos,
             c, d] for a, b, c, d in self._cs
        ]

    def get_contig(self):
        return self._contig

    def modify_contig(self, contig):
        old_contig = self._contig
        self._contig = contig
        return old_contig

    def modify_start_pos(self, start_pos):
        old_start_pos = self._start_pos
        self._start_pos = start_pos
        return old_start_pos

    def get_start_pos(self):
        return self._start_pos

    def get_strand(self):
        return self._strand

    def _get_read_length(self):
        read_len = 0
        for a, b, c, d in self.get_relative_position():
            if c == ":":
                read_len += int(d)
            elif c == "~":
                pass
            elif c == "+":
                read_len += len(d)
            elif c == "-":
                pass
            elif c == "*":
                read_len += 1
        return read_len

    def get_read_length(self):
        return self._read_length

    def _get_element_read_location(self):
        read_location = list()
        now_len = 0
        previous_len = 0
        for a, b, c, d in self.get_relative_position():
            if c == ":":
                now_len += int(d)
                read_location.append([
                    previous_len, now_len
                ])
            elif c == "~":
                read_location.append([
                    previous_len, now_len
                ])
            elif c == "+":
                now_len += len(d)
                read_location.append([
                    previous_len, now_len
                ])
            elif c == "-":
                read_location.append([
                    previous_len, now_len
                ])
            elif c == "*":
                now_len += 1
                read_location.append([
                    previous_len, now_len
                ])
            else:
                pass
            previous_len = now_len
        return read_location

    def get_read_location(self):
        strand = self.get_strand()
        read_loc = self._get_element_read_location()
        relative_pos = self.get_relative_position()
        read_len = self.get_read_length()
        new_cs = list()
        for a, b in zip(read_loc, relative_pos):
            if strand == '+':
                new_cs.append(a + b[2:])
            else:
                new_cs.append(
                    [read_len - a[1],
                     read_len - a[0]] + b[2:]
                )
        return new_cs

    def get_normalized_read_location(self):
        read_loc = self.get_read_location()
        read_len = self.get_read_length()
        normalized_loc = [
            [a / read_len, b / read_len, c, d]
            for a, b, c, d in read_loc
        ]
        return normalized_loc

    def _count_element_in_read_location_bin(self, mark):
        assert mark in ['~', '*', ':', '+', '-'], \
            "Mark should be in [~, *, :, +, -]."
        normalized_loc = self.get_normalized_read_location()
        bin_counter = Counter()
        for loc in normalized_loc:
            if loc[2] != mark:
                continue
            else:
                if loc[0] >= 0 and loc[0] < 0.25:
                    bin_counter['[0, 0.25)'] += 1
                elif loc[0] >= 0.25 and loc[0] < 0.5:
                    bin_counter['[0.25, 0.5)'] += 1
                elif loc[0] >= 0.5 and loc[0] < 0.75:
                    bin_counter['[0.5, 0.75)'] += 1
                else:
                    bin_counter['[0.75, 1]'] += 1
        return bin_counter

    def count_insertion_in_read_location_bin(self):
        return self._count_element_in_read_location_bin('+')

    def count_deletion_in_read_location_bin(self):
        return self._count_element_in_read_location_bin('-')

    def count_mismatch_in_read_location_bin(self):
        return self._count_element_in_read_location_bin('*')

    def count_intron_in_read_location_bin(self):
        return self._count_element_in_read_location_bin('~')

    def _get_marks(self, mark, coordinate):
        assert mark in ['~', '*', ':', '+', '-'], \
            "Mark should be in [~, *, :, +, -]."
        assert coordinate in [
            'contig', 'relative_contig', 'read', 'normalized_read'
        ], 'coordinate should be in [contig, relative_contig, read, normalized_read]'
        if coordinate == 'contig':
            outlist = [
                a for a in self.get_contig_position()
                if a[2] == mark
            ]
        elif coordinate == 'relative_contig':
            outlist = [
                a for a in self.get_relative_position()
                if a[2] == mark
            ]
        elif coordinate == 'read':
            outlist = [
                a for a in self.get_read_location()
                if a[2] == mark
            ]
        elif coordinate == 'normalized_read':
            outlist = [
                a for a in self.get_normalized_read_location()
                if a[2] == mark
            ]
        else:
            pass
        return outlist

    def get_introns(self, coordinate = 'contig'):
        assert coordinate in [
            'contig', 'relative_contig', 'read', 'normalized_read'
        ], 'coordinate should be in [contig, relative_contig, read, normalized_read]'
        return self._get_marks('~', coordinate = coordinate)

    def get_mismatches(self, coordinate = 'contig'):
        assert coordinate in [
            'contig', 'relative_contig', 'read', 'normalized_read'
        ], 'coordinate should be in [contig, relative_contig, read, normalized_read]'
        return self._get_marks('*', coordinate = coordinate)

    def get_matches(self, coordinate = 'contig'):
        assert coordinate in [
            'contig', 'relative_contig', 'read', 'normalized_read'
        ], 'coordinate should be in [contig, relative_contig, read, normalized_read]'
        return self._get_marks(':', coordinate = coordinate)

    def get_insertions(self, coordinate = 'contig'):
        assert coordinate in [
            'contig', 'relative_contig', 'read', 'normalized_read'
        ], 'coordinate should be in [contig, relative_contig, read, normalized_read]'
        return self._get_marks('+', coordinate = coordinate)

    def get_deletions(self, coordinate = 'contig'):
        assert coordinate in [
            'contig', 'relative_contig', 'read', 'normalized_read'
        ], 'coordinate should be in [contig, relative_contig, read, normalized_read]'
        return self._get_marks('-', coordinate = coordinate)

    def get_intron_count(self):
        return self._intron_count

    def _count_intron(self):
        return len(self.get_introns())

    def get_splice_pair_count(self):
        return self._splice_pair_count

    def _count_splice_pair(self):
        introns = self.get_introns()
        splice_pairs = [
            re.sub('[0-9]+', '-', a[3])
            for a in introns
        ]
        pair_ctr = Counter()
        for a in splice_pairs:
            pair_ctr[a] += 1
        return pair_ctr

    def get_splice_site_count(self):
        return (self._splice_l_count,
                self._splice_r_count)

    def _count_splice_site(self):
        introns = self.get_introns()
        splice_sites = [
            re.sub('[0-9]+', ' ', a[3]).strip().split()
            for a in introns
        ]
        l_ctr = Counter()
        r_ctr = Counter()
        for a, b in splice_sites:
            l_ctr[a] += 1
            r_ctr[b] += 1
        return l_ctr, r_ctr

    def get_mismatch_count(self):
        return self._mismatch_count

    def _count_mismatch(self):
        mismatches = self.get_mismatches()
        return len(mismatches)

    def get_mismatch_type_count(self):
        return self._mismatch_type_count

    def _count_mismatch_type(self):
        mis = self.get_mismatches()
        mis_ctr = Counter()
        for a in [item[3] for item in mis]:
            mis_ctr['{}{}'.format(a[0], a[1])] += 1
        return mis_ctr

    def get_insertion_count(self):
        return self._insertion_count

    def get_insertion_length(self):
        return self._insertion_length

    def _count_insertion(self):
        insertion = self.get_insertions()
        return len(insertion)

    def _calculate_insertion_length(self):
        insertion = self.get_insertions()
        if len(insertion) > 0:
            return sum(
                [len(a[3]) for a in insertion]
            )
        else:
            return 0

    def get_deletion_count(self):
        return self._deletion_count

    def get_deletion_length(self):
        return self._deletion_length

    def _count_deletion(self):
        deletion = self.get_deletions()
        return len(deletion)

    def _calculate_deletion_length(self):
        deletion = self.get_deletions()
        if len(deletion) > 0:
            return sum(
                [len(a[3]) for a in deletion]
            )
        else:
            return 0

    def __str__(self):
        outstring = 'CS tag: {} introns, {} mismatches, {} indels'.format(
            self._intron_count,
            sum(self._mismatch_count.values()),
            self._insertion_count + self._deletion_count
        )
        return outstring

    def __repr__(self):
        outstring = '\n'.join([
            "CS tag:",
            "  {} introns".format(self._intron_count),
            "  {} mismatches".format(sum(self._mismatch_count)),
            "  {} indels".format(self._insertion_count + self._deletion_count),
            "--------------------",
            self.get_cs_tag_string()
        ])
        return outstring

########################################
