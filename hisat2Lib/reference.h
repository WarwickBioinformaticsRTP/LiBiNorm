/*
 * Copyright 2011, Ben Langmead <langmea@cs.jhu.edu>
 *
 * This file is part of Bowtie 2.
 *
 * Bowtie 2 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Bowtie 2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Bowtie 2.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef REFERENCE_H_
#define REFERENCE_H_

#include <vector>
#include <string>
#include <stdexcept>
#include <fcntl.h>
#include <sys/stat.h>
#include <utility>
#ifdef BOWTIE_MM
#include <sys/mman.h>
#include <sys/shm.h>
#endif
#include "ref_read.h"
//#include "sstring.h"
#include "btypes.h"


/**
 * Concrete reference representation that bulk-loads the reference from
 * the bit-pair-compacted binary file and stores it in memory also in
 * bit-pair-compacted format.  The user may request reference
 * characters either on a per-character bases or by "stretch" using
 * getBase(...) and getStretch(...) respectively.
 *
 * Most of the complexity in this class is due to the fact that we want
 * to represent references with ambiguous (non-A/C/G/T) characters but
 * we don't want to use more than two bits per base.  This means we
 * need a way to encode the ambiguous stretches of the reference in a
 * way that is external to the bitpair sequence.  To accomplish this,
 * we use the RefRecords vector, which is stored in the .3.ebwt index
 * file.  The bitpairs themselves are stored in the .4.ebwt index file.
 *
 * Once it has been loaded, a BitPairReference is read-only, and is
 * safe for many threads to access at once.
 */
class BitPairReference {

public:
	/**
	 * Load from .3.ebwt/.4.ebwt Bowtie index files.
	 */
	BitPairReference(
		const string& in,
		bool sanity = false,
		bool verbose = false,
		bool startVerbose = false);

	~BitPairReference();

	/**
	 * Load a stretch of the reference string into memory at 'dest'.
	 *
	 * This implementation scans linearly through the records for the
	 * unambiguous stretches of the target reference sequence.  When
	 * there are many records, binary search would be more appropriate.
	 */
	int getStretchNaive(
		uint32_t *destU32,
		size_t tidx,
		size_t toff,
		size_t count) const;

	/**
	 * Load a stretch of the reference string into memory at 'dest'.
	 *
	 * This implementation scans linearly through the records for the
	 * unambiguous stretches of the target reference sequence.  When
	 * there are many records, binary search would be more appropriate.
	 */
	int getStretch(
		uint32_t *destU32,
		size_t tidx,
		size_t toff,
		size_t count
		ASSERT_ONLY(, vector<uint32_t>& destU32_2)) const;

	/**
	 * Return the number of reference sequences.
	 */
	TIndexOffU numRefs() const {
		return nrefs_;
	}
protected:

	uint32_t byteToU32_[256];

	vector<RefRecord> recs_;       /// records describing unambiguous stretches
	// following two lists are purely for the binary search in getStretch
	vector<TIndexOffU> cumUnambig_; // # unambig ref chars up to each record
	vector<TIndexOffU> cumRefOff_;  // # ref chars up to each record
	vector<TIndexOffU> refLens_;    /// approx lens of ref seqs (excludes trailing ambig chars)
	vector<TIndexOffU> refOffs_;    /// buf_ begin offsets per ref seq
	vector<TIndexOffU> refRecOffs_; /// record begin/end offsets per ref seq
	uint8_t *buf_;      /// the whole reference as a big bitpacked byte array
	TIndexOffU bufSz_;    /// size of buf_
	TIndexOffU bufAllocSz_;
	TIndexOffU nrefs_;    /// the number of reference sequences
	bool     loaded_;   /// whether it's loaded
	bool     sanity_;   /// do sanity checking
	bool     verbose_;
	ASSERT_ONLY(vector<uint32_t> tmp_destU32_);
};

#endif
