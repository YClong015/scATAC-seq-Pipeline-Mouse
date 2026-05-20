#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --job-name=DWLD_CNX0739875
#SBATCH --time=48:00:00
#SBATCH --partition=general
#SBATCH --account=a_nefzger
#SBATCH --output=DWLD_%j.out
#SBATCH --error=DWLD_%j.err

# --- 变量定义 ---

# 1. 你要递归下载的 FTP 目录 URL
FTP_DIR_URL="ftp://ftp2.cngb.org/pub/CNSA/data6/CNP0004399/CNS0809369/CNX0739880/CNR0841136/"
# 2. wget 将在 $TMPDIR 中创建的目录名
DOWNLOAD_DIR_NAME="CNX0739880"

# 3. 最终数据存放的目标目录
DATADIR=/QRISdata/Q8448/Mouse_disease_data/Lung

# --- 脚本执行 ---

echo "--- 开始任务：下载 $FTP_DIR_URL ---"
echo "--- 临时目录: $TMPDIR"
echo "--- 目标目录: $DATADIR"

# 切换到 HPC 节点的本地临时存储 ($TMPDIR)
cd $TMPDIR

echo "--- 正在启动 wget 递归下载..."

# 执行 wget 递归下载
# 这将在 $TMPDIR 中创建 CNX0739875 目录
wget -c -nH -np -r -R "index.html*" --cut-dirs 5 "$FTP_DIR_URL"

echo "--- wget 下载完成 ---"

# 检查下载目录是否真的被创建了
if [ ! -d "$DOWNLOAD_DIR_NAME" ]; then
    echo "错误：wget 未能成功创建目录 $DOWNLOAD_DIR_NAME"
    exit 1
fi

echo "--- 下载完成。正在将整个目录复制到 $DATADIR..."

# 递归复制 (cp -r) 整个下载目录到你的永久存储
# 这将在 DATADIR 下创建一个名为 CNX0739875 的新目录
cp -r ${DOWNLOAD_DIR_NAME} ${DATADIR}/

echo "--- 任务全部完成。数据位于: ${DATADIR}/${DOWNLOAD_DIR_NAME} ---"
