import os
import os.path
import datetime
import operator

import Bio.SeqIO

samples_dates = []

align_fname = '/idi/sabeti-scratch/swohl/mumps/outbreaker/v3-final/II-outbreak-v3.aligned.pruned.fasta'

for seq_rec in Bio.SeqIO.parse(align_fname, 'fasta'):
    fields = seq_rec.id.split('|')
    sample = fields[0]
    samples_dates.append((sample, datetime.datetime.strptime(fields[-1], '%Y-%m-%d')))
    Bio.SeqIO.write(seq_rec, 'data/02_assembly/{}.fasta'.format(sample), 'fasta-2line')

print(len(samples_dates))
earliest_date = min(map(operator.itemgetter(1), samples_dates))
print(earliest_date)
for sample, date in samples_dates:
    fname = 'tmp/06_outbreak/{}.days_since_outbreak_start.txt'.format(sample)
    days = (date - earliest_date).days
    print(fname, days)
    with open(fname, 'w') as out:
        out.write(str(days))
    assert os.path.isfile('tmp/cleaned/{}.raw.cleaned.bam'.format(sample)), 'missing for: {}'.format(sample)
    os.symlink(os.path.abspath(os.path.realpath('tmp/cleaned/{}.raw.cleaned.bam'.format(sample))), 
               'data/01_per_sample/{}.cleaned.bam'.format(sample))

with open('samples-assembly.txt', 'w') as out:
    out.write('\n'.join(map(operator.itemgetter(0), samples_dates)))

