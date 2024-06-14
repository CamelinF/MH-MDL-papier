/*
 * The MIT License (MIT)
 *
 * Copyright (c) 2014 University of Bonn
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package de.unibonn.realkd.patterns.models.table;

import java.util.Comparator;
import java.util.List;

import com.google.common.collect.ImmutableList;

import ua.ac.be.mime.errorhandling.NotYetImplementedException;

public class BinaryOrdinalDimension<T> implements Dimension {

	private static final ImmutableList<String> BIN_CAPTIONS = ImmutableList.of("low", "high");
	
	private final Comparator<T> comparator;
	private final T threshold;

	BinaryOrdinalDimension(Comparator<T> comparator, T threshold) {
		this.comparator = comparator;
		this.threshold = threshold;
	}

	@Override
	public int bin(Object value) {
		if (comparator.compare((T) value, threshold)<=0) {
			return 0;
		}
		return 1;
	}

	@Override
	public List<String> binCaptions() {
		return BIN_CAPTIONS;
	}

	@Override
	public int numberOfBins() {
		return 2;
	}

	@Override
	public Dimension merge(int to, int from) {
		return null;
	}

	@Override
	public Dimension condition(int[] binsToCondition) {
		return null;
	}
	
}