function plotECGSegmented(states, tmAsynch, sigAsynch, seqViterbi, marker)


plot(tmAsynch, sigAsynch);
hold on
state = 0;
for ll = 1 : length(seqViterbi)
    if states == 4
        [state1Idx, ~] = find(seqViterbi{1, ll} == 1);
        [state2Idx, ~] = find(seqViterbi{1, ll} == 2);
        [state3Idx, ~] = find(seqViterbi{1, ll} == 3);
        [state4Idx, ~] = find(seqViterbi{1, ll} == 4);
        state1Idx = state + state1Idx;
        state2Idx = state + state2Idx;
        state3Idx = state + state3Idx;
        state4Idx = state + state4Idx;
        plot(tmAsynch(state1Idx), sigAsynch(state1Idx), marker{1})
        plot(tmAsynch(state2Idx), sigAsynch(state2Idx), marker{2})
        plot(tmAsynch(state3Idx), sigAsynch(state3Idx), marker{3})
        plot(tmAsynch(state4Idx), sigAsynch(state4Idx), marker{4})
    elseif states == 6
        [state1Idx, ~] = find(seqViterbi{1, ll} == 1);
        [state2Idx, ~] = find(seqViterbi{1, ll} == 2);
        [state3Idx, ~] = find(seqViterbi{1, ll} == 3);
        [state4Idx, ~] = find(seqViterbi{1, ll} == 4);
        [state5Idx, ~] = find(seqViterbi{1, ll} == 5);
        [state6Idx, ~] = find(seqViterbi{1, ll} == 6);
        state1Idx = state + state1Idx;
        state2Idx = state + state2Idx;
        state3Idx = state + state3Idx;
        state4Idx = state + state4Idx;
        state5Idx = state + state5Idx;
        state6Idx = state + state6Idx;
        plot(tmAsynch(state1Idx), sigAsynch(state1Idx), marker{1})
        plot(tmAsynch(state2Idx), sigAsynch(state2Idx), marker{2})
        plot(tmAsynch(state3Idx), sigAsynch(state3Idx), marker{3})
        plot(tmAsynch(state4Idx), sigAsynch(state4Idx), marker{4})
        plot(tmAsynch(state5Idx), sigAsynch(state5Idx), marker{5})
        plot(tmAsynch(state6Idx), sigAsynch(state6Idx), marker{6})
    elseif states == 3
        [state1Idx, ~] = find(seqViterbi{1, ll} == 1);
        [state2Idx, ~] = find(seqViterbi{1, ll} == 2);
        [state3Idx, ~] = find(seqViterbi{1, ll} == 3);
        state1Idx = state + state1Idx;
        state2Idx = state + state2Idx;
        state3Idx = state + state3Idx;
        plot(tmAsynch(state1Idx), sigAsynch(state1Idx), marker{1})
        plot(tmAsynch(state2Idx), sigAsynch(state2Idx), marker{2})
        plot(tmAsynch(state3Idx), sigAsynch(state3Idx), marker{3})
    end
    
    state = length(seqViterbi{1, ll}) + state;
    
end
hold off
xlabel('tmAsynch')
ylabel('sigAsynch')

if states == 4
    legend('signal','idle', 'P-wave', 'QRS-complex', 'T-wave')
elseif states == 6
    legend('signal', 'idle1', 'P-wave', 'idle2', 'QRS-complex', 'idle3', 'T-wave')
elseif states == 3
    legend('signal', 'P-wave', 'QRS-complex', 'T-wave')
end
end