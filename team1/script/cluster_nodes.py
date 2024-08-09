import networkx as nx
import sys

# 读取文件并创建图
G = nx.read_edgelist(sys.argv[1])

# 计算连通分量
clusters = list(nx.connected_components(G))

# 创建一个字典存储节点到cluster的映射
node_to_cluster = {}
for cluster_id, cluster in enumerate(clusters):
    for node in cluster:
        node_to_cluster[node] = cluster_id

# 输出结果
with open(sys.argv[2], 'w') as f:
    for node, cluster_id in sorted(node_to_cluster.items()):
        f.write(f"{node} {cluster_id}\n")

