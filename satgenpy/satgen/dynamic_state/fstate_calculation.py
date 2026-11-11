import math
import networkx as nx # FOR floyd_warshall_numpy

# 1. Calculate shortest path without the use of GS as nodes
def calculate_fstate_shortest_path_without_gs_relaying(
        output_dynamic_state_dir,
        time_since_epoch_ns,
        num_satellites,
        num_ground_stations,
        sat_net_graph_only_satellites_with_isls,
        num_isls_per_sat,
        gid_to_sat_gsl_if_idx,
        ground_station_satellites_in_range_candidates,
        sat_neighbor_to_if,
        prev_fstate,
        enable_verbose_logs
):

    # Calculate shortest path distances
    if enable_verbose_logs:
        print("  > Calculating Floyd-Warshall for graph without ground-station relays")
    # (Note: Numpy has a deprecation warning here because of how networkx uses matrices)
    dist_sat_net_without_gs = nx.floyd_warshall_numpy(sat_net_graph_only_satellites_with_isls) # minimum distance between satellite a and b passing only through ISL

    # Forwarding state
    fstate = {}

    # Now write state to file for complete graph
    output_filename = output_dynamic_state_dir + "/fstate_" + str(time_since_epoch_ns) + ".txt"
    if enable_verbose_logs:
        print("  > Writing forwarding state to: " + output_filename)
    with open(output_filename, "w+") as f_out:

        # SATELLITES TO GROUND STATION
        # From the satellites attached to the destination ground station, select the one which promises the shortest path to the destination ground station (getting there + last hop)
        dist_satellite_to_ground_station = {}
        for curr in range(num_satellites):
            for dst_gid in range(num_ground_stations):
                dst_gs_node_id = num_satellites + dst_gid # identify final GS node (satellites  i=0...n_sat-1 )
                                                          #                        (ground stations i=n_sat...n_sat+n_gs-1) 

                # Among the satellites in range of the destination ground station,
                # find the one which promises the shortest distance
                possible_dst_sats = ground_station_satellites_in_range_candidates[dst_gid] # satellites that can see the gs
                possibilities = []
                # for each satellite b(1), calculates distance ISL from curr to b(1). Adds final GSL distance
                for b in possible_dst_sats:
                    if not math.isinf(dist_sat_net_without_gs[(curr, b[1])]):  # Must be reachable
                        possibilities.append(
                            (
                                dist_sat_net_without_gs[(curr, b[1])] + b[0],
                                b[1]
                            )
                        )
                possibilities = list(sorted(possibilities)) # contains pairs of info possibilities(i) = (total_distance, dst_sat). it is sorted and the best possibilkity is chosen

                # By default, if there is no satellite in range for the
                # destination ground station, it will be dropped (indicated by -1)
                next_hop_decision = (-1, -1, -1)
                distance_to_ground_station_m = float("inf")
                if len(possibilities) > 0: # if there is a valid possibility, choose 1st on the list (it has been sorted). 1st satellite offers the lowest cost to GS
                    dst_sat = possibilities[0][1]
                    distance_to_ground_station_m = possibilities[0][0]

                    # If the current node is not the satellite of best possibility, determine how to get to the satellite
                    if curr != dst_sat:

                        # Among its neighbors, find the one which promises the lowest distance to reach the destination satellite (iterate through neighbours, keep the option if the distance is minimized)
                        best_distance_m = 1000000000000000
                        for neighbor_id in sat_net_graph_only_satellites_with_isls.neighbors(curr):
                            distance_m = (
                                    sat_net_graph_only_satellites_with_isls.edges[(curr, neighbor_id)]["weight"]
                                    +
                                    dist_sat_net_without_gs[(neighbor_id, dst_sat)]
                            )
                            if distance_m < best_distance_m:
                                next_hop_decision = (
                                    neighbor_id,
                                    sat_neighbor_to_if[(curr, neighbor_id)],
                                    sat_neighbor_to_if[(neighbor_id, curr)]
                                )
                                best_distance_m = distance_m

                    # If current node is the final satellite
                    else:
                        # This is the destination satellite, as such the next hop is the ground station itself
                        next_hop_decision = (
                            dst_gs_node_id,
                            num_isls_per_sat[dst_sat] + gid_to_sat_gsl_if_idx[dst_gid],
                            0
                        )

                # In any case, save the distance of the satellite to the ground station to re-use when we calculate ground station to ground station forwarding
                dist_satellite_to_ground_station[(curr, dst_gs_node_id)] = distance_to_ground_station_m

                # Write to forwarding state (only if next_hop decision changed)
                if not prev_fstate or prev_fstate[(curr, dst_gs_node_id)] != next_hop_decision:
                    f_out.write("%d,%d,%d,%d,%d\n" % (
                        curr,
                        dst_gs_node_id,
                        next_hop_decision[0],
                        next_hop_decision[1],
                        next_hop_decision[2]
                    ))
                fstate[(curr, dst_gs_node_id)] = next_hop_decision

        # GROUND STATIONS TO GROUND STATIONS (if I'm in the initial GS (source), how do I get to the final one (destination))
        # Choose the source satellite which promises the shortest path
        for src_gid in range(num_ground_stations):
            for dst_gid in range(num_ground_stations):
                if src_gid != dst_gid:
                    # IDs of GS nodes
                    src_gs_node_id = num_satellites + src_gid
                    dst_gs_node_id = num_satellites + dst_gid

                    # Among the satellites in range of the source ground station, find the one which promises the shortest distance
                    possible_src_sats = ground_station_satellites_in_range_candidates[src_gid]
                    possibilities = []
                    for a in possible_src_sats:
                        best_distance_offered_m = dist_satellite_to_ground_station[(a[1], dst_gs_node_id)]
                        if not math.isinf(best_distance_offered_m):
                            possibilities.append(
                                (
                                    a[0] + best_distance_offered_m,
                                    a[1]
                                )
                            )
                    possibilities = sorted(possibilities)

                    # By default, if there is no satellite in range for one of the ground stations, it will be dropped (indicated by -1)
                    next_hop_decision = (-1, -1, -1)
                    if len(possibilities) > 0:
                        src_sat_id = possibilities[0][1]
                        next_hop_decision = (
                            src_sat_id,
                            0,
                            num_isls_per_sat[src_sat_id] + gid_to_sat_gsl_if_idx[src_gid]
                        )

                    # Update forwarding state
                    if not prev_fstate or prev_fstate[(src_gs_node_id, dst_gs_node_id)] != next_hop_decision:
                        f_out.write("%d,%d,%d,%d,%d\n" % (
                            src_gs_node_id,
                            dst_gs_node_id,
                            next_hop_decision[0],
                            next_hop_decision[1],
                            next_hop_decision[2]
                        ))
                    fstate[(src_gs_node_id, dst_gs_node_id)] = next_hop_decision

    # Finally return result
    return fstate

# 2. Calculate shortest path with the use of GS (GS can be an intermediate node of the communication route)
def calculate_fstate_shortest_path_with_gs_relaying(
        output_dynamic_state_dir,
        time_since_epoch_ns,
        num_satellites,
        num_ground_stations,
        sat_net_graph,
        num_isls_per_sat,
        gid_to_sat_gsl_if_idx,
        sat_neighbor_to_if,
        prev_fstate,
        enable_verbose_logs
):

    # Calculate shortest paths
    if enable_verbose_logs:
        print("  > Calculating Floyd-Warshall for graph including ground-station relays")
    # (Note: Numpy has a deprecation warning here because of how networkx uses matrices)
    dist_sat_net = nx.floyd_warshall_numpy(sat_net_graph) # sat_net_graph includes  all nodes and possible links

    # Forwarding state
    fstate = {}

    # Now write state to file for complete graph
    output_filename = output_dynamic_state_dir + "/fstate_" + str(time_since_epoch_ns) + ".txt"
    if enable_verbose_logs:
        print("  > Writing forwarding state to: " + output_filename)
    with open(output_filename, "w+") as f_out:

        # Satellites and ground stations to ground stations
        for current_node_id in range(num_satellites + num_ground_stations):
            for dst_gid in range(num_ground_stations):
                dst_gs_node_id = num_satellites + dst_gid # id final GS

                # Cannot forward to itself
                if current_node_id != dst_gs_node_id:

                    # Among its neighbors, find the one which promises the lowest distance to reach the destination satellite
                    next_hop_decision = (-1, -1, -1)
                    best_distance_m = 1000000000000000
                    for neighbor_id in sat_net_graph.neighbors(current_node_id):

                        # Any neighbor must be reachable
                        if math.isinf(dist_sat_net[(current_node_id, neighbor_id)]):
                            raise ValueError("Neighbor cannot be unreachable")

                        # Calculate distance = next-hop + distance the next hop node promises
                        distance_m = (
                            sat_net_graph.edges[(current_node_id, neighbor_id)]["weight"]
                            +
                            dist_sat_net[(neighbor_id, dst_gs_node_id)]
                        )
                        if (
                                not math.isinf(dist_sat_net[(neighbor_id, dst_gs_node_id)])
                                and
                                distance_m < best_distance_m
                        ):

                            # Check node identifiers to determine what are the correct interface identifiers
                            if current_node_id >= num_satellites and neighbor_id < num_satellites:  # GS to sat.
                                my_if = 0 # GS uses interface 0
                                next_hop_if = ( # Satellite uses interface GSL that corresponds to that GS
                                    num_isls_per_sat[neighbor_id]
                                    +
                                    gid_to_sat_gsl_if_idx[current_node_id - num_satellites]
                                )

                            elif current_node_id < num_satellites and neighbor_id >= num_satellites:  # Sat. to GS
                                my_if = ( # Satellite uses interface GSL that corresponds to that GS
                                    num_isls_per_sat[current_node_id]
                                    +
                                    gid_to_sat_gsl_if_idx[neighbor_id - num_satellites]
                                )
                                next_hop_if = 0  # GS uses interface 0

                            elif current_node_id < num_satellites and neighbor_id < num_satellites:  # Sat. to sat.
                                my_if = sat_neighbor_to_if[(current_node_id, neighbor_id)] # Satellite uses mapped interface 
                                next_hop_if = sat_neighbor_to_if[(neighbor_id, current_node_id)] # Satellite uses mapped interface 

                            else:  # GS to GS
                                raise ValueError("GS-to-GS link cannot exist")

                            # Write the next-hop decision
                            next_hop_decision = (
                                neighbor_id,  # Next-hop node identifier
                                my_if,        # My outgoing interface id
                                next_hop_if   # Next-hop incoming interface id
                            )

                            # Update best distance found
                            best_distance_m = distance_m

                    # Write to forwarding state
                    if not prev_fstate or prev_fstate[(current_node_id, dst_gs_node_id)] != next_hop_decision:
                        f_out.write("%d,%d,%d,%d,%d\n" % (
                            current_node_id,
                            dst_gs_node_id,
                            next_hop_decision[0],
                            next_hop_decision[1],
                            next_hop_decision[2]
                        ))
                    fstate[(current_node_id, dst_gs_node_id)] = next_hop_decision

    # Finally return result
    return fstate
