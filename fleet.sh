#!/bin/bash
# Generic fleet controller: runs a plan produced by `pi <x> --plan [alpha opts]`.
#
#   pi 1e24 --plan --alpha-fit=A,B > plan24.txt
#   AMI=... SG=... KEY=... PEM=... bash fleet.sh plan24.txt ~/fleet24
#
# The plan carries the parameter contract (X/Y/SEGW) and the task specs
# (kind.a:b/N); this script is pure transport: launch tasks under a vCPU cap,
# collect fragments over ssh, terminate agents, merge with --check.
# State = fragment files + instance tags — kill and rerun at any time.
set -u
export PATH="$HOME/.local/bin:$PATH"

PLAN=${1:?usage: fleet.sh <planfile> <workdir>}
WORK=${2:?usage: fleet.sh <planfile> <workdir>}
AMI=${AMI:-ami-0fd23c7420775ce6c}
SG=${SG:-sg-00229d4999c53c8b3}
KEY=${KEY:-prime-count}
TYPE=${TYPE:-c8g.xlarge}
CAP=${CAP:-7}
PEM=${PEM:-$HOME/.ec2/prime-count.pem}
HASH=${HASH:-$(cd "$(dirname "$0")" && git rev-parse HEAD)}
MERGE_T=${MERGE_T:-6}

X=$(awk '/^X /{print $2}' "$PLAN")
Y=$(awk '/^Y /{print $2}' "$PLAN")
SEGW=$(awk '/^SEGW /{print $2}' "$PLAN")
mapfile -t SPECS < <(grep -E '^[AB]\.' "$PLAN")
[ -n "$X" ] && [ -n "$Y" ] && [ -n "$SEGW" ] && [ "${#SPECS[@]}" -gt 0 ] || {
  echo "bad plan file" >&2; exit 2; }
RUNTAG="fleet-$X-${HASH:0:7}"
mkdir -p "$WORK"
cd "$WORK"
echo "plan: x=$X y=$Y segw=$SEGW tasks=${#SPECS[@]} code=$HASH tag=$RUNTAG"

fragfile() { local s="${1//[:\/.]/_}"; echo "frag_${s}.out"; }

launch() { # $1 = kind.a:b/N
  local kind="${1%%.*}" iv="${1#*.}" flag="--blocks"
  [ "$kind" = "A" ] && flag="--asig"
  cat > ud.sh <<EOF
#!/bin/bash
sudo -u ubuntu bash -c '
set -e
cd /home/ubuntu/prime-counting
git fetch origin
git checkout $HASH
/home/ubuntu/.local/bin/zig build-exe -O ReleaseFast -mcpu=native --dep common --dep rs -Mroot=src/combinatorial/pi.zig -Mcommon=src/common.zig -Mrs=src/rangesieve.zig -femit-bin=./pi-task
./pi-task $X --y $Y --segw $SEGW $flag $iv --emit /home/ubuntu/frag.out -t 4 -v
touch /home/ubuntu/DONE
'
EOF
  aws ec2 run-instances --image-id "$AMI" --instance-type "$TYPE" \
    --key-name "$KEY" --security-group-ids "$SG" \
    --instance-initiated-shutdown-behavior terminate \
    --user-data file://ud.sh \
    --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=$RUNTAG-$kind-${iv//[:\/]/-}},{Key=pi-x,Value=$RUNTAG},{Key=pi-blocks,Value=$1},{Key=pi-code,Value=$HASH}]" \
    --query 'Instances[0].InstanceId' --output text 2>&1
}

NTOT=${#SPECS[@]}
while :; do
  roster=$(aws ec2 describe-instances \
    --filters "Name=tag:pi-x,Values=$RUNTAG" "Name=instance-state-name,Values=pending,running" \
    --query 'Reservations[].Instances[].[InstanceId,Tags[?Key==`pi-blocks`]|[0].Value,PublicIpAddress]' \
    --output text)
  nrun=0
  live=""
  while read -r id spec ip; do
    [ -z "${id:-}" ] && continue
    nrun=$((nrun + 1))
    live="$live $spec"
    f=$(fragfile "$spec")
    [ -f "$f" ] && continue
    if [ "$ip" != "None" ] && ssh -n -i "$PEM" -o StrictHostKeyChecking=no -o ConnectTimeout=5 \
        "ubuntu@$ip" test -f /home/ubuntu/DONE 2>/dev/null; then
      scp -q -i "$PEM" -o StrictHostKeyChecking=no "ubuntu@$ip:/home/ubuntu/frag.out" "$f" &&
        echo "$(date +%H:%M:%S) collected $spec from $id" &&
        aws ec2 terminate-instances --instance-ids "$id" --query 'TerminatingInstances[0].CurrentState.Name' --output text >/dev/null &&
        nrun=$((nrun - 1))
    fi
  done <<< "$roster"

  ndone=0
  next=""
  for sp in "${SPECS[@]}"; do
    if [ -f "$(fragfile "$sp")" ]; then
      ndone=$((ndone + 1))
    elif [ -z "$next" ] && ! grep -qF " $sp" <<< "$live "; then
      next="$sp"
    fi
  done
  echo "$(date +%H:%M:%S) status: $ndone/$NTOT collected, $nrun running"
  [ "$ndone" -eq "$NTOT" ] && break

  while [ -n "$next" ] && [ "$nrun" -lt "$CAP" ]; do
    out=$(launch "$next")
    case "$out" in
      i-*) echo "$(date +%H:%M:%S) launched $out <- $next" ;;
      *) echo "$(date +%H:%M:%S) launch deferred ($next): ${out:0:80}"; break ;;
    esac
    live="$live $next"
    nrun=$((nrun + 1))
    next=""
    for sp in "${SPECS[@]}"; do
      [ -f "$(fragfile "$sp")" ] && continue
      if ! grep -qF " $sp" <<< "$live "; then next="$sp"; break; fi
    done
  done
  sleep 30
done

echo "=== all $NTOT fragments — merge ==="
"$(dirname "$0")/pi" --merge frag_*.out --check -t "$MERGE_T" --pin
echo "=== FLEET RUN COMPLETE ==="
