import numpy as np
from functools import reduce

def rotate_points_forpiece(points, indices, rotate):
    #input : N x 2 
    indices = np.argmax(indices,1)
    indices[indices==0] = 1000
    unique_indices = np.unique(indices)
    num_unique_indices = len(unique_indices)
    rotated_points = np.zeros_like(points)
    rotation_angles = []
    for i in unique_indices:
        idx = (indices == i)
        selected_points = points[idx] 
        #조각의 point 1개/  k x 2
        #rotation_degree = 0 if i==1 else (np.random.rand() * 360)
        # rotation_angle = 0 
        # rotation_angle = 0 if i==0 else (np.random.randint(4) * 90)
        rotation_angle = rotate
        rotation_matrix = np.array([
            [np.cos(rotation_angle), -np.sin(rotation_angle)], # this is selected for return
            [np.sin(rotation_angle), np.cos(rotation_angle)]
        ])

        # 반시계 방향 회전 
        rotated_selected_points = np.matmul(rotation_matrix, selected_points.T).T
        rotated_points[idx] = rotated_selected_points
        # rotation_matrix[0,1] = 1 if rotation_angle<np.pi else -1
                                # 1 x 2  ->  k x 2
        rotation_angles.extend(rotation_matrix[0:1].repeat(rotated_selected_points.shape[0], axis=0))
    return rotated_points, rotation_angles

def find_canonicalize(poly):
    import numpy as np
    vectors = poly
    # Cartesian 좌표를 Polar 좌표로 변환
    polar_coordinates = np.array([[np.sqrt(x**2 + y**2), np.arctan2(y, x)] for x, y in vectors])
    for rho, theta in polar_coordinates:
        print(f"Rho: {rho:.2f}, Theta: {theta:.2f} radians")

    # 벡터들의 좌표
    vectors = poly

    # 벡터들의 각도를 구하고, 이를 x-y 좌표계에 매핑
    x = np.cos(np.arctan2(vectors[:, 1], vectors[:, 0]))
    y = np.sin(np.arctan2(vectors[:, 1], vectors[:, 0]))

    # 평균 각도 계산
    mean_x = np.mean(x)
    mean_y = np.mean(y)

    if mean_x == 0 and mean_y == 0:
        print("The mean vector is at the origin. The mean angle is undefined.")
    else:
        mean_angle = np.arctan2(mean_y, mean_x)
        print(f"Mean Angle: {mean_angle:.2f} radians")
    
    canonical_theta=mean_angle
    return polar_coordinates, canonical_theta

def raddif(rada, radb):
    import math
    res=rada-radb
    if (res > math.pi):
        return (res - 2*math.pi)
    elif (res< -1*math.pi):
        return (res + 2*math.pi)
    else:
        return (res)

class ordered_cyclic():
    def __init__(self, nump):
        self.nump=nump[nump[:,1].argsort()]
        
        self.len=nump.shape[0]
    def distance(self, indices):
        xx=[]
        for i in range(len(indices)-1):
            xx.append([indices[i+1]-indices[i], raddif(self.nump[indices[i+1]][1],self.nump[indices[i]][1])])
        xx.append([self.len-indices[-1]+indices[0], raddif(self.nump[indices[0]][1],self.nump[indices[-1]][1])])

        return np.array(xx) 
    def indexsame(self, indexa, indexb, threshold=1e-6):
        return (((indexa-indexb)**2).mean(0)<threshold)
    def displacementsame(self, r1, r2, threshold=1e-8):
        return (((r1-r2)**2).mean(0)<threshold)
    def group(self):
         # 2차원 좌표 (x, y)의 배열
        coordinates = self.nump

        # 허용 오차
        tolerance = 1e-6

        # x 좌표를 그룹화
        groups = {}
        for z, coord in enumerate(coordinates):
            x, y = coord
            matched = False
            for key in groups.keys():
                if np.isclose(x, key, atol=tolerance):
                    groups[key].append(z)
                    matched = True
                    break
            
            if not matched:
                groups[x] = [z]

        #key(거리에 따라 오름차순 정렬)
        groups = dict(sorted(groups.items()))
        #dict를 오름차순 정렬.
        # 결과 출력
        for key, values in groups.items():
            print(f"x ~ {key}: {values}")
        
        NewGroup={}
        for key in groups.keys():
            A=groups[key]
            Te=self.distance(A)
            NewGroup[key]=Te
        def Findcycle(coordinates ,threshold=1e-5):
            import numpy as np
            from scipy.spatial.distance import cdist
            # 모든 좌표 쌍 간의 유클리디안 거리 계산
            distances = cdist(coordinates, coordinates)

            # 각 좌표와의 거리가 임계값 이내인 좌표 찾기
            similar_points = distances < threshold
            # 반복성 검사
            iscycle=False
            for i in range(1, len(coordinates)//2+1):
                if len(coordinates)%i!=0:
                    continue
                
                pattern = similar_points[:i]
                zz=True
                for t, pat in enumerate(pattern):
                    l=t
                    while l < len(pat):
                        if pat[l]==False:
                            zz=False
                        l+=i
                
                if zz:
                    iscycle=i
                    print(f"Pattern repeats every {i} coordinates.")
                    break
            if iscycle==False:
                return find_parsed(coordinates)
            else: 
                return None
        
        for key in NewGroup.keys():
             #거리가 같은것들을 모아서 차이를 기반으로 (순서 유지)된 list에서 cycle을 찾는다. 
             #+ 해당 coordinates를 이용해 검색, deterministic한 결과를 뽑아낸다
             A=Findcycle(NewGroup[key])
             if A!=None:
                 print(groups[key][A])
                 return NewGroup, groups, groups[key][A]
             else:
                 continue
        for key in NewGroup.keys():
             return NewGroup, groups, groups[key][0]
    def process(self):
        Groups=self.group()
        #순서 유지

def find_parsed(original_array):
    print(original_array)
    original_array=np.array(original_array)
    # 두 좌표 간의 거리를 계산
    from scipy.spatial.distance import cdist
    distances = cdist(original_array, original_array)
    # 거리 임계값을 설정 (이 값을 조절하여 클러스터의 크기를 변경할 수 있음)
    threshold = 1e-6
    # 클러스터 리스트 초기화
    clusters = []
    clusters_id=[]
    # 각 좌표에 대해
    for i in range(len(original_array)):
        # 해당 좌표가 이미 클러스터에 포함되어 있는지 확인
        if any(i in cluster for cluster in clusters):
            continue
        
        # 해당 좌표와 다른 모든 좌표 간의 거리를 확인하여 임계값보다 작은 거리를 가진 좌표들의 인덱스 찾기
        near_points = np.where(distances[i] < threshold)[0]
        
        #------------near point 삭제해야됨 

        # 클러스터에 추가
        clusters_id.append(i)
        clusters.append(list(near_points))
    
    '''print(original_array)
    print(clusters)
    '''
    #p=sorted(list(set(np.ndarray.tolist(original_array))))
    original_order=list([list(original_array[i]) for i in clusters_id])
    p= sorted(list([list(original_array[i]) for i in clusters_id]))
    matching=[p.index(original) for original in original_order]
    print("-------")
    print(p)
    print(original_order)
    print(matching)
    new_clusters_id=[]
    for k in range(len(matching)):
        new_clusters_id.append(clusters_id[matching[k]])
    #clusters_id=[clusters_id[matching[k]] for k in range(len(matching))]
    new_clusters=[]
    for k in range(len(matching)):
        new_clusters.append(clusters[matching.index(k)])
    #clusters_id=[clusters_id[matching[k]] for k in range(len(matching))]
    print(clusters)
    # clusters & clusters_id 순서 변경해야됨 
    NewAssign=np.zeros_like(original_array[:, 0])
    for k, cluster in enumerate(new_clusters):
        for index in cluster:
            NewAssign[index]=k
    
    print(new_clusters)

    for k, item in enumerate(p): # id(1, 0.2324) 으로 검색
        print(item)
        K=find_unique_parsed_arrays(NewAssign, item, new_clusters[k])
        if K!=None:
            print(K)
            return K[1]
    return None
        
def find_unique_parsed_arrays(original_array, delimiter, cluster):
    # 특정 원소를 기준으로 배열 파싱
    #delimiter_indices = np.where(original_array == delimiter)[0]
    delimiter_indices=np.array(cluster)
    parse=[]
    for k in range(len(delimiter_indices)-1):
        parse.append(np.concatenate((original_array[delimiter_indices[k]:delimiter_indices[k+1]], np.array([delimiter_indices[k]])), axis=0) )
    parse.append(np.concatenate((original_array[delimiter_indices[-1]:], original_array[:delimiter_indices[0]], np.array([delimiter_indices[-1]])), axis=0))
    parsed_arrays = parse
    
    # 파싱한 배열들을 오름차순으로 정렬
    tempray=[parsed_arrays[i][:-1] for i in range(len(parsed_arrays))]
    sums = [np.sum(arr) for arr in tempray]
    # 합을 기준으로 배열들을 정렬
    sorted_indices = np.argsort(sums, axis=0)
 
    # 정렬된 배열 출력
    sorted_arrays = [parsed_arrays[i][:-1] for i in sorted_indices]
    sorted_finalindice=[parsed_arrays[i][-1] for i in sorted_indices]
    # 결과 출력
    '''for arr in sorted_arrays:
        print(arr)'''
    '''import pdb
    pdb.set_trace()'''
    sorted_parsed_arrays = sorted_arrays
    print("------------sorted------------")
    print(sorted_arrays)
    for ll, parsei in enumerate(sorted_parsed_arrays):
        cnt=0
        for p in sorted_parsed_arrays:
            if np.ndarray.tolist(parsei)==np.ndarray.tolist(p):
                cnt+=1
        if cnt==1:
            subarray_length = len(parsei)

            # 특정 배열의 시작 인덱스를 찾기
            start_index = -1
            return parsei, int(sorted_finalindice[ll])
    return None
    # Unique한 배열 찾기
    unique_arrays = np.unique(sorted_parsed_arrays, axis=0)
    import pdb
    pdb.set_trace()
    # Unique한 배열에 대해 해당 배열을 파싱하는 데에 사용된 특정 원소의 원본 배열에서의 인덱스 찾기
    for unique_array in unique_arrays:
        for i, parsed_array in enumerate(sorted_parsed_arrays):
            if np.array_equal(unique_array, parsed_array):
                print(f"Unique array: {unique_array}, Index of delimiter in original array: {delimiter_indices[i]}")
                break

                
        
        
                


    
    #input: rho_theta_coordinate( N x 2 ) -> rho-theta orientation
def find_equivariant_canonical(polar):
    A=ordered_cyclic(rho_theta_coordinates)
    R, U , targetarrow = A.group()
    return targetarrow


if __name__=="__main__": 
    A=ordered_cyclic(rho_theta_coordinates)
    R, U , targetarrow = A.group()