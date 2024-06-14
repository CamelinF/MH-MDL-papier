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

package de.unibonn.realkd.data.table.attribute;

import static java.util.stream.Collectors.toList;

import java.util.Date;
import java.util.List;

import com.fasterxml.jackson.annotation.JsonCreator;
import com.fasterxml.jackson.annotation.JsonProperty;

import de.unibonn.realkd.common.base.Identifier;

/**
 * 
 * @author Sandy Moens
 * 
 * @since 0.3.0
 * 
 * @version 0.4.0
 * 
 */
public class DefaultDateAttribute extends DefaultAttribute<Date> {

//	private final CategoryContainer<Date> categoryContainer;
	@JsonProperty("values")
	private List<Long> values;

	@JsonCreator
	DefaultDateAttribute(@JsonProperty("identifier") Identifier identifier, @JsonProperty("name") String name,
			@JsonProperty("description") String description, @JsonProperty("values") List<Date> values) {
		super(identifier, name, description, values, Date.class);
		this.values = values.stream().map(d -> d.getTime()).collect(toList());
//		this.categoryContainer = new CategoryContainer<>(nonMissingValues());
	}

//	@Override
//	public List<Double> categoryFrequencies() {
//		return categoryContainer.getCategoryFrequencies();
//	}
//
//	@Override
//	public List<Date> categories() {
//		return categoryContainer.getCategories();
//	}

}
