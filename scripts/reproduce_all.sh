#!/bin/bash
# reproduce_all.sh — One-command reproduction of all BrM results
# Seed: 42, matches manuscript. Requires Modal credentials.
set -e

echo "=== evo2-e2e Reproduction Script ==="
echo "Timestamp: $(date -u +%Y%m%dT%H%M%SZ)"
echo "Seed: 42"
echo ""

# Step 1: Stress tests (verify anti-hallucination)
echo "── Step 1: Anti-hallucination stress tests ──"
python tests/stress/test_antihallucination.py --fast
echo ""

# Step 2: Brain met pipeline (fast mode)
echo "── Step 2: BrM Target-Lock pipeline (fast) ──"
python brain_met/pipeline/run_brm_pipeline.py --fast --seed 42
echo ""

# Step 3: BBB Transit step specifically
echo "── Step 3: BBB Transit scoring ──"
python brain_met/pipeline/run_brm_pipeline.py --step bbb_transit --seed 42
echo ""

echo "=== Reproduction complete. Results in data/brain_met/ ==="
