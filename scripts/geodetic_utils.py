from obspy.geodetics.base import gps2dist_azimuth
from obspy.geodetics import locations2degrees


def add_distance_backazimuth_to_df(df, event):
    backazimuth_list = []
    distance_list = []
    for _, row in df.iterrows():
        distance, azimuth, backazimuth = gps2dist_azimuth(
            lat1=row.latitude,
            lon1=row.longitude,
            lat2=event["latitude"],
            lon2=event["longitude"],
        )
        distance = locations2degrees(
            lat1=row.latitude,
            long1=row.longitude,
            lat2=event["latitude"],
            long2=event["longitude"],
        )
        # Append results to lists
        backazimuth_list.append(backazimuth)
        distance_list.append(distance)
    # Add the new columns to the DataFrame
    df["backazimuth"] = backazimuth_list
    df["distance"] = distance_list
    return df
