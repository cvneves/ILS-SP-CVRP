#ifndef INFOSEQ
#define INFOSEQ

#include "Data.h"

/** Stores information about a subsequence of a path, which 
 * is used to Evaluate movements efficiently */
struct Subsequence
{
    double duration;
    double timewarp;
    double earliest;
    double latest;
    double distance;
    double load;
#ifdef ROBUST_CONCATENATION_TEST
    std::vector<double> critLoad;
    std::vector<double> critTimewarp;
    std::vector<double> critEarliest;
    std::vector<double> critLatest;
    std::vector<double> critDuration;
    int nb_clients;
#endif
    double cap_violation;
    int first_cl; // index of the first client in the subsequence
    int last_cl;  // index of the last client in the subsequence

    inline static Subsequence Concatenate(Data *data, Subsequence &subseq1, Subsequence &subseq2)
    {
        Subsequence subseq;

        double delta, deltaWT, deltaTW;
        delta = subseq1.duration - subseq1.timewarp + data->time_cost[subseq1.last_cl][subseq2.first_cl];
        deltaWT = std::max(subseq2.earliest - delta - subseq1.latest, 0.0);
        deltaTW = std::max(subseq1.earliest + delta - subseq2.latest, 0.0);
        subseq.duration = subseq1.duration + subseq2.duration + data->time_cost[subseq1.last_cl][subseq2.first_cl] + deltaWT;
        subseq.timewarp = subseq1.timewarp + subseq2.timewarp + deltaTW;
        subseq.earliest = std::max(subseq2.earliest - delta, subseq1.earliest) - deltaWT;
        subseq.latest = std::min(subseq2.latest - delta, subseq1.latest) + deltaTW;
        subseq.distance = subseq1.distance + subseq2.distance + data->time_cost[subseq1.last_cl][subseq2.first_cl];
        subseq.load = subseq1.load + subseq2.load;
        subseq.first_cl = subseq1.first_cl;
        subseq.last_cl = subseq2.last_cl;
        subseq.cap_violation = std::max(subseq.load - data->capacity, 0.0);

#ifdef ROBUST_CONCATENATION_TEST
        subseq.nb_clients = subseq1.nb_clients + subseq2.nb_clients;
        subseq.critLoad = std::vector<double>(data->budget_demand + 1);
        subseq.critLoad[0] = 0;
        for (int gamma = 1, i1 = 0, i2 = 0; gamma <= data->budget_demand; gamma++)
        {
            if (subseq1.critLoad[i1 + 1] + subseq2.critLoad[i2] > subseq1.critLoad[i1] + subseq2.critLoad[i2 + 1])
                i1++;
            else
                i2++;
            subseq.critLoad[gamma] = subseq1.critLoad[i1] + subseq2.critLoad[i2];
        }

        subseq.critTimewarp = std::vector<double>(data->budget_time + 1);
        subseq.critTimewarp[0] = subseq.timewarp;
        subseq.critEarliest = std::vector<double>(data->budget_time + 1);
        subseq.critEarliest[0] = subseq.earliest;
        subseq.critLatest = std::vector<double>(data->budget_time + 1);
        subseq.critLatest[0] = subseq.latest;
        subseq.critDuration = std::vector<double>(data->budget_time + 1);
        subseq.critDuration[0] = subseq.duration;

        double critTime[2] = {data->time_cost[subseq1.last_cl][subseq2.first_cl],
                              data->time_cost[subseq1.last_cl][subseq2.first_cl] +
                                  data->time_deviation[subseq1.last_cl][subseq2.first_cl]};
        double critDelta, critDeltaWT, critDeltaTW;

        for (int gamma = 1, i1 = 0, i2 = 0, i3 = 0; gamma <= data->budget_time; gamma++)
        {
            i1++;
            bool inc1 = true, inc2 = false, inc3 = false;
            critDelta = subseq1.critDuration[i1] - subseq1.critTimewarp[i1] + critTime[i3];
            critDeltaWT = std::max(subseq2.critEarliest[i2] - critDelta - subseq1.critLatest[i1], 0.0);
            critDeltaTW = std::max(subseq1.critEarliest[i1] + critDelta - subseq2.critLatest[i2], 0.0);
            subseq.critDuration[gamma] = subseq1.critDuration[i1] + subseq2.critDuration[i2] + critTime[i3] + critDeltaWT;
            subseq.critTimewarp[gamma] = subseq1.critTimewarp[i1] + subseq2.critTimewarp[i2] + critDeltaTW;
            subseq.critEarliest[gamma] = std::max(subseq2.critEarliest[i2] - critDelta, subseq1.critEarliest[i1]) - critDeltaWT;
            subseq.critLatest[gamma] = std::min(subseq2.critLatest[i2] - critDelta, subseq1.critLatest[i1]) + critDeltaTW;
            i1--;

            i2++;
            critDelta = subseq1.critDuration[i1] - subseq1.critTimewarp[i1] + critTime[i3];
            critDeltaWT = std::max(subseq2.critEarliest[i2] - critDelta - subseq1.critLatest[i1], 0.0);
            critDeltaTW = std::max(subseq1.critEarliest[i1] + critDelta - subseq2.critLatest[i2], 0.0);
            double currDuration = subseq1.critDuration[i1] + subseq2.critDuration[i2] + critTime[i3] + critDeltaWT;
            double currTimewarp = subseq1.critTimewarp[i1] + subseq2.critTimewarp[i2] + critDeltaTW;
            double currEarliest = std::max(subseq2.critEarliest[i2] - critDelta, subseq1.critEarliest[i1]) - critDeltaWT;
            double currLatest = std::min(subseq2.critLatest[i2] - critDelta, subseq1.critLatest[i1]) + critDeltaTW;

            if (currDuration > subseq.critDuration[gamma])
            {
                subseq.critDuration[gamma] = currDuration;
                subseq.critTimewarp[gamma] = currTimewarp;
                subseq.critEarliest[gamma] = currEarliest;
                subseq.critLatest[gamma] = currLatest;
                inc1 = false, inc2 = true, inc3 = false;
            }
            i2--;

            if (i3 == 1)
                continue;

            i3++;

            critDelta = subseq1.critDuration[i1] - subseq1.critTimewarp[i1] + critTime[i3];
            critDeltaWT = std::max(subseq2.critEarliest[i2] - critDelta - subseq1.critLatest[i1], 0.0);
            critDeltaTW = std::max(subseq1.critEarliest[i1] + critDelta - subseq2.critLatest[i2], 0.0);
            currDuration = subseq1.critDuration[i1] + subseq2.critDuration[i2] + critTime[i3] + critDeltaWT;
            currTimewarp = subseq1.critTimewarp[i1] + subseq2.critTimewarp[i2] + critDeltaTW;
            currEarliest = std::max(subseq2.critEarliest[i2] - critDelta, subseq1.critEarliest[i1]) - critDeltaWT;
            currLatest = std::min(subseq2.critLatest[i2] - critDelta, subseq1.critLatest[i1]) + critDeltaTW;
            if (currDuration > subseq.critDuration[gamma])
            {
                subseq.critDuration[gamma] = currDuration;
                subseq.critTimewarp[gamma] = currTimewarp;
                subseq.critEarliest[gamma] = currEarliest;
                subseq.critLatest[gamma] = currLatest;
                inc1 = false, inc2 = false, inc3 = true;
            }
            i3--;

            i1 += inc1;
            i2 += inc2;
            i3 += inc3;
        }
#endif

        return subseq;

        // if (subseq2.earliest - delta > subseq1.latest)
        // {
        //     deltaTW = 0;
        //     subseq.earliest = subseq1.latest;
        //     subseq.latest = subseq1.latest;
        //     deltaWT = subseq2.earliest - delta - subseq1.latest;
        // }
        // else if (subseq1.earliest + delta > subseq2.latest)
        // {
        //     deltaWT = 0;
        //     subseq.earliest = subseq1.earliest;
        //     subseq.latest = subseq1.earliest;
        //     deltaTW = subseq1.earliest + delta - subseq2.latest;
        // }
        // else
        // {
        //     deltaWT = 0;
        //     deltaTW = 0;
        //     subseq.earliest = std::max(subseq2.earliest - delta, subseq1.earliest);
        //     subseq.latest = std::min(subseq2.latest - delta, subseq1.latest);
        // }

        // subseq.duration = subseq1.duration + subseq2.duration + data->time_cost[subseq1.last_cl][subseq2.first_cl] + deltaWT;
        // subseq.timewarp = subseq1.timewarp + subseq2.timewarp + deltaTW;
        // subseq.distance = subseq1.distance + subseq2.distance + data->time_cost[subseq1.last_cl][subseq2.first_cl];
        // subseq.load = subseq1.load + subseq2.load;
        // subseq.first_cl = subseq1.first_cl;
        // subseq.last_cl = subseq2.last_cl;
        // subseq.cap_violation = std::max(subseq.load - data->capacity, 0.0);
    }
};

#endif
