---
# SPDX-FileCopyrightText: 2017-2023 CERN and the Corryvreckan authors
# SPDX-License-Identifier: CC-BY-4.0 OR MIT
---
# FilterEvents
**Maintainer**: Lennart Huth (lennart.huth@desy.de)  
**Module Type**: *GLOBAL*  
**Status**: Work in progress  

### Description
This module allows for a filtering of events based on the number of tracks/clusters on the clipboard. If no value for the parameters is given, the filter on the parameter is inactive.
Events can be skipped based on the number of tracks and clusters per plane. The accepted cluster range is defined for all planes, that are not a `DUT` or a `Auxiliary`.
Events can also be skipped based on tag values where the tag value is either checked to be witihin a specified range, or matching at least one value of a specified list.

### Parameters

* `exclude_trigger_windows`: List of trigger windows to exclude in analysis. A trigger window is a list of a lower and a higher trigger ID bound. For example, the trigger window `[51, 99]` would ignore all events with trigger numbers between 50 and 100. Defaults to an empty list.
* `min_tracks`: Minimum number of tracks to continue analyzing event. Per default inactive.
* `max_tracks`: Maximum number of tracks to continue analyzing event. Per default inactive.
* `only_tracks_on_dut`: Apply the `min_tracks` and `max_tracks` cut only on tracks that intersect the DUT. If more than one DUT is defined in the geometry,the bool is set to its default. Defaults to `false`
* `min_clusters_per_plane`: Minimum number of clusters on each plane to continue analyzing event. Planes of type `DUT` and `auxiliary` are excluded in this cut. Per default inactive.
* `max_clusters_per_plane`: Maximum number of clusters on each plane to continue analyzing event. Planes of type `DUT` and `auxiliary` are excluded in this cut. Per default inactive..
* `filter_tags` : List of tag value requirements to be fulfilled to continue analyzing event. The list of tag-filter pairs should be enclosed in global pair of brackets `[ list_of_tag_filter_pairs ]`. Each individual tag-filter pair should be enclosed in local pair of square brackets `[]`, and each element enclosed in a pair of quotes `""` (ie : `["tag_name", "tag_filter"]`). If the tag filter is to be a range of values (casted to `double` type) the minimum and maximum values should be enclosed in a pair of squared brackets `[]` and separated by a semicolomn `:` (ie : `["tag_name", "[value_min:value_max]"]`). If the tag filter is to be a list of values, the values should separated by a comma `,` (ie : `["tag_name", "value1,value2,...,valueN"]`). For range-based filtering, the values are compared as doubles and events with tag values outside the range are not processed. For list-based filtering, the values are compared as strings and events with tag values different from the listed ones are not processed.
* `min_event_duration`: Minimal event length to continue analyzing the event. Per default inactive.
* `max_event_duration`: Maximal event length to continue analyzing the event. Per default inactive.

### Plots produced
* 1D histogram of the events filtered by a certain parameter and the total number of filtered events.

### Usage
```toml
[FilterEvents]
exclude_trigger_windows = [[0, 100], [2600, 2600]]
min_tracks = 0
max_tracks = 10
min_clusters_per_plane = 0
max_clusters_per_plane = 100
min_event_duration = 0ns
max_event_duration = 50ns
```
