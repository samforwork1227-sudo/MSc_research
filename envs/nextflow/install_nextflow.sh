#!/usr/bin/env bash
set -euo pipefail

PARENT_DIR="/scratch/users/k25053715/MSc_research/0_project_env"
BASE_DIR="${PARENT_DIR}/nextflow_env"

INSTALL_DIR="${BASE_DIR}/bin"
NXF_HOME_DIR="${BASE_DIR}/nxf_home"
APPTAINER_CACHE_DIR="${BASE_DIR}/apptainer_cache"
APPTAINER_TMP_DIR="${BASE_DIR}/apptainer_tmp"
JDK_DIR="${BASE_DIR}/jdk"

mkdir -p "${BASE_DIR}"
mkdir -p "${INSTALL_DIR}" "${NXF_HOME_DIR}" "${APPTAINER_CACHE_DIR}" "${APPTAINER_TMP_DIR}" "${JDK_DIR}"

echo "==> Nextflow environment directory:"
echo "    ${BASE_DIR}"
echo "==> Java will be installed here:"
echo "    ${JDK_DIR}"

JDK_TARBALL="${BASE_DIR}/OpenJDK17U-jdk_x64_linux_hotspot_17.0.16_8.tar.gz"

if [ ! -f "${JDK_TARBALL}" ]; then
    echo "==> Downloading JDK 17 (183MB, please wait...)"
    wget -O "${JDK_TARBALL}" \
      "https://github.com/adoptium/temurin17-binaries/releases/download/jdk-17.0.16%2B8/OpenJDK17U-jdk_x64_linux_hotspot_17.0.16_8.tar.gz"
fi

echo "==> Found JDK tarball:"
echo "    ${JDK_TARBALL}"

echo "==> Cleaning old JDK directory"
rm -rf "${JDK_DIR}"
mkdir -p "${JDK_DIR}"

echo "==> Unpacking JDK 17"
tar -xzf "${JDK_TARBALL}" -C "${JDK_DIR}" --strip-components=1

export PATH="${JDK_DIR}/bin:${PATH}"
export JAVA_HOME="${JDK_DIR}"
export NXF_JAVA_HOME="${JDK_DIR}"

echo "==> Checking Java"
java -version

# Download Nextflow
echo "==> Downloading Nextflow"
cd "${BASE_DIR}"
curl -fsSL https://get.nextflow.io | bash

chmod +x nextflow
mv -f nextflow "${INSTALL_DIR}/nextflow"

cat > "${BASE_DIR}/load_nextflow_env.sh" <<EOF
#!/usr/bin/env bash
export PATH="${JDK_DIR}/bin:${INSTALL_DIR}:\$PATH"
export JAVA_HOME="${JDK_DIR}"
export NXF_JAVA_HOME="${JDK_DIR}"
export NXF_HOME="${NXF_HOME_DIR}"
export APPTAINER_CACHEDIR="${APPTAINER_CACHE_DIR}"
export APPTAINER_TMPDIR="${APPTAINER_TMP_DIR}"
EOF

chmod +x "${BASE_DIR}/load_nextflow_env.sh"

export PATH="${INSTALL_DIR}:${PATH}"
export NXF_HOME="${NXF_HOME_DIR}"
export APPTAINER_CACHEDIR="${APPTAINER_CACHE_DIR}"
export APPTAINER_TMPDIR="${APPTAINER_TMP_DIR}"

echo "==> Verifying Nextflow"
"${INSTALL_DIR}/nextflow" -version

echo ""
echo "==> Installation complete"
echo "To use this environment, run:"
echo "  source ${BASE_DIR}/load_nextflow_env.sh"
echo ""
echo "Then test:"
echo "  which java"
echo "  java -version"
echo "  which nextflow"
echo "  nextflow -version"
echo "  echo \$NXF_HOME"
